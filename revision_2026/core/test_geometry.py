"""Analytical geometry and conservation checks for the scientific remapper."""
import unittest
import numpy as np
from affine import Affine
from common import overlaps,exact_joint_support,separable_resample

class GeometryTests(unittest.TestCase):
    def test_extensive_split_and_partial_coverage(self):
        src=np.array([0.,1.,2.]);dst=np.array([-.5,.5,1.5,2.5])
        out=overlaps(src,dst,'source')@np.array([10.,30.])
        np.testing.assert_allclose(out,[5,20,15]);self.assertAlmostEqual(out.sum(),40)
        out=overlaps(src,dst,'target')@np.ones(2)
        np.testing.assert_allclose(out,[.5,1,.5])
    def test_joint_masks_not_product_of_means(self):
        tr=Affine(1,0,0,0,-1,1)
        cc=np.array([[255,0]],dtype='uint8');cv=np.ones_like(cc,dtype=bool)
        uv=np.array([[0,1]],dtype=bool)
        out=exact_joint_support(cc,cv,tr,uv,tr,np.array([0.,2.]),np.array([0.,1.]),[.3])
        self.assertAlmostEqual(float(out['unw'][0,0]),.5)
        self.assertAlmostEqual(float(out['q30'][0,0]),0.)
    def test_half_pixel_shift_intersection(self):
        cc=np.array([[255,0]],dtype='uint8');cv=np.ones_like(cc,dtype=bool)
        uv=np.array([[1,0]],dtype=bool)
        out=exact_joint_support(cc,cv,Affine(1,0,0,0,-1,1),uv,Affine(1,0,.5,0,-1,1),np.array([0.,2.]),np.array([0.,1.]),[.3])
        self.assertAlmostEqual(float(out['q30'][0,0]),.25)
    def test_outside_coverage_and_valid_zero_coherence(self):
        a=np.array([[0]],dtype='uint8');v=np.ones_like(a,dtype=bool)
        out=exact_joint_support(a,v,Affine(1,0,0,0,-1,1),v,Affine(1,0,0,0,-1,1),np.array([-1.,0.,1.,2.]),np.array([0.,1.]),[0.,.3])
        np.testing.assert_allclose(out['q00'],[[0,1,0]])
        np.testing.assert_allclose(out['q30'],[[0,0,0]])
    def test_spherical_mass_transfer(self):
        a=np.array([[1.,2.],[3.,4.]])
        b=separable_resample(a,[0,1,2],[10,11,12],[0,.3,.8,1.4,2],[10,10.2,11.2,12],normalization='source')
        self.assertAlmostEqual(b.sum(),10.,places=12)

if __name__=='__main__':unittest.main()
