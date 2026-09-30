from pathlib import Path
import os,json,hashlib
import numpy as np
from scipy.sparse import csr_matrix

CODE_ROOT=Path(__file__).resolve().parent
ROOT=Path(os.environ.get('INSAR_CORE_ROOT', CODE_ROOT)).resolve()
CFG=json.loads((CODE_ROOT/'config.json').read_text(encoding='utf-8'))

def resolve_data_root(config_key, env_key, default):
    value=os.environ.get(env_key, CFG.get(config_key, default))
    path=Path(value).expanduser()
    return (path if path.is_absolute() else ROOT/path).resolve()

# Resolve source-data roots relative to this checkout by default. Environment
# variables let users point to separately retained legacy/input archives.
PROJECT=resolve_data_root('project_root','INSAR_PROJECT_ROOT','..')
BACKUP=resolve_data_root('backup_root','INSAR_BACKUP_ROOT','../backup_inputs')
OUT=ROOT/'results';OUT.mkdir(parents=True,exist_ok=True)

def resolve_recorded_path(value):
    """Locate a retained manifest input after a checkout/dataset move.

    Historical manifest strings remain unchanged as provenance. Only known
    project/backup prefixes are rebased; missing or ambiguous inputs fail.
    Callers must still retain/verify the recorded content checksum.
    """
    recorded=Path(value)
    if recorded.exists():return recorded
    text=str(value).replace('\\','/')
    candidates=[]
    if '/strengthening/' in text:
        candidates.append(ROOT/text.split('/strengthening/',1)[1])
    if text.startswith('D:/InSAR_revision_20260926/'):
        candidates.append(PROJECT/text[len('D:/InSAR_revision_20260926/'):])
    if '/nature/' in text:
        candidates.append(BACKUP/text.split('/nature/',1)[1])
    if not recorded.is_absolute():candidates.append(ROOT/recorded)
    existing={p.resolve() for p in candidates if p.is_file()}
    if len(existing)!=1:
        raise FileNotFoundError(f'Cannot uniquely resolve recorded input {value!r}; '
            'configure INSAR_PROJECT_ROOT and INSAR_BACKUP_ROOT and verify the input manifest.')
    return existing.pop()

def dump(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False),encoding='utf-8')

def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    return h.hexdigest()

def overlaps(source_edges,target_edges,normalization='target'):
    """Sparse interval intersections. Edges strictly increasing.

    For latitude pass sin(latitude), for longitude pass radians. Target
    normalization averages an intensive field over the WHOLE target interval,
    including uncovered parts as zero. Source normalization transfers extensive
    mass conservatively. No extrapolation or data-dependent renormalization.
    """
    source_edges=np.asarray(source_edges,dtype='float64');target_edges=np.asarray(target_edges,dtype='float64')
    if np.any(np.diff(source_edges)<=0) or np.any(np.diff(target_edges)<=0):raise ValueError('Nonincreasing edges')
    rows=[];cols=[];values=[]
    for i,(a,b) in enumerate(zip(target_edges[:-1],target_edges[1:])):
        j=max(0,int(np.searchsorted(source_edges,a,side='right')-1))
        while j<len(source_edges)-1 and source_edges[j]<b:
            width=min(b,source_edges[j+1])-max(a,source_edges[j])
            if width>0:
                denom=(b-a) if normalization=='target' else (source_edges[j+1]-source_edges[j]) if normalization=='source' else 1.
                rows.append(i);cols.append(j);values.append(width/denom)
            j+=1
    return csr_matrix((values,(rows,cols)),shape=(len(target_edges)-1,len(source_edges)-1))

def separable_resample(a,source_x,source_y,target_x,target_y,normalization='target',spherical=True):
    """All y edges and array rows ordered south-to-north (or ascending northing)."""
    if spherical:
        source_x=np.deg2rad(source_x);target_x=np.deg2rad(target_x)
        source_y=np.sin(np.deg2rad(source_y));target_y=np.sin(np.deg2rad(target_y))
    wx=overlaps(source_x,target_x,normalization);wy=overlaps(source_y,target_y,normalization)
    return (wx @ (wy @ a).T).T

def edges(transform,shape):
    h,w=shape
    if abs(transform.b)>1e-12 or abs(transform.d)>1e-12 or transform.a<=0 or transform.e>=0:raise ValueError('Unsupported raster transform')
    return transform.c+np.arange(w+1)*transform.a, (transform.f+np.arange(h+1)*transform.e)[::-1]

def pixel_areas(x,y):
    return (6371008.8**2)*np.diff(np.sin(np.deg2rad(y)))[:,None]*np.diff(np.deg2rad(x))[None,:]

def merged_edges(a,b,lo,hi):
    x=np.sort(np.concatenate([[lo,hi],a[(a>lo)&(a<hi)],b[(b>lo)&(b<hi)]]))
    return x[np.r_[True,np.diff(x)>1e-11]]

def sample_array(array,transform,x,y,invalid=0):
    cc=np.floor((x-transform.c)/transform.a).astype(int)
    rr=np.floor((y-transform.f)/transform.e).astype(int)
    inside_y=(rr>=0)&(rr<array.shape[0]);inside_x=(cc>=0)&(cc<array.shape[1])
    sampled=array[np.ix_(np.clip(rr,0,array.shape[0]-1),np.clip(cc,0,array.shape[1]-1))].copy()
    sampled[~inside_y,:]=invalid;sampled[:,~inside_x]=invalid
    return sampled

def exact_joint_support(cc,cc_valid,cc_tr,unw_valid,unw_tr,tx,ty,thresholds):
    """Exact spherical rectangular overlap of two north-up geographic masks.

    Split at BOTH source grids' edges before intersection; no multiplication
    of separately averaged fractional masks. Output rows are south-to-north.
    """
    cx,cy=edges(cc_tr,cc.shape);ux,uy=edges(unw_tr,unw_valid.shape)
    sx=merged_edges(cx,ux,tx[0],tx[-1]);sy=merged_edges(cy,uy,ty[0],ty[-1])
    xx=(sx[:-1]+sx[1:])/2;yy=(sy[:-1]+sy[1:])/2
    c=sample_array(cc,cc_tr,xx,yy)
    cvalid=sample_array(cc_valid,cc_tr,xx,yy).astype(bool)
    u=sample_array(unw_valid,unw_tr,xx,yy).astype(bool)
    wx=overlaps(np.deg2rad(sx),np.deg2rad(tx),'target')
    wy=overlaps(np.sin(np.deg2rad(sy)),np.sin(np.deg2rad(ty)),'target')
    def transfer(mask):
        dst=(wx @ (wy @ mask.astype('float32')).T).T
        if dst.min() < -1e-8 or dst.max()>1+1e-8:raise ValueError('Invalid area fraction')
        return np.clip(dst,0,1).astype('float32')
    result={'unw':transfer(u)}
    for q in thresholds:
        result[f'q{int(q*100):02d}']=transfer(u & cvalid & (c>=int(np.ceil(q*255-1e-10))))
    return result
