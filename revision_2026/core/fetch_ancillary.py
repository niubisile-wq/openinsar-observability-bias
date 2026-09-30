"""Retrieve fixed ancillary inputs and retain their exact local checksums."""
from common import ROOT,CFG,dump
from fetch_external import fetch

def main():
    records=[]
    sources=[('https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_N12E099_Map.tif',ROOT/'external/worldcover2021_N12E099.tif')]
    base=CFG['catalogue_url'].replace('interferograms/','metadata/')
    for name in [CFG['frame']+'.geo.hgt.tif',CFG['frame']+'.geo.U.tif','baselines','metadata.txt']:
        sources.append((base+name,ROOT/'external/licsar_metadata'/name))
    for url,path in sources:
        path.parent.mkdir(parents=True,exist_ok=True);record=fetch(url,path);record['relative_path']=str(path.relative_to(ROOT)).replace('\\','/');records.append(record)
    dump(ROOT/'external/ancillary_manifest.json',records)
    print('ANCILLARY INPUTS RECORDED',len(records),flush=True)

if __name__=='__main__':main()
