"""Regression tests for fully missing and zero-weight observed units."""
from pathlib import Path
import sys,unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'paired_support'))
from masking_estimators import geometry,retained_statistics,evaluate

class AvailabilityTests(unittest.TestCase):
    def calculate(self,pop,missing):
        pop=np.array([pop],float);rate=np.array([[-10,0,-10,0]],float)
        valid=np.ones_like(pop,bool);missing=np.array([missing],bool)
        g=geometry(pop,valid,2)
        return evaluate(g,rate,missing,-5,retained_statistics(g,rate,missing,True),True)

    def test_entire_missing_unit_is_undefined(self):
        for row in self.calculate([90,10,30,70],[True,True,False,False]):
            self.assertEqual(row['unestimated_allocation_units'],100)
            self.assertEqual(row['common_computable_allocation_units'],100)
            self.assertEqual(row['common_signed_error_pp'],0)
            self.assertEqual(row['whole_domain_lower_bound_pct'],15)
            self.assertEqual(row['whole_domain_upper_bound_pct'],65)

    def test_zero_population_observation_does_not_define_prevalence(self):
        rows=self.calculate([10,0,30,70],[True,False,False,False])
        self.assertEqual(rows[0]['unestimated_allocation_units'],10)
        self.assertEqual(rows[1]['unestimated_allocation_units'],0)
        for row in rows:
            self.assertEqual(row['common_computable_allocation_units'],100)
            self.assertEqual(row['common_signed_error_pp'],0)

    def test_no_computable_domain_has_no_fabricated_error(self):
        for row in self.calculate([90,10,30,70],[True]*4):
            self.assertEqual(row['unestimated_allocation_share_pct'],100)
            self.assertTrue(np.isnan(row['common_signed_error_pp']))
            self.assertTrue(np.isnan(row['common_relative_error_pct']))
            self.assertEqual(row['whole_domain_bound_width_pp'],100)

    def test_defined_but_unrepresentative_observations(self):
        for row in self.calculate([90,10,30,70],[True,False,True,False]):
            self.assertEqual(row['common_signed_error_pp'],-60)
            self.assertEqual(row['common_local_L1_pp'],60)
            self.assertEqual(row['whole_domain_upper_bound_pct'],60)

if __name__=='__main__':unittest.main()
