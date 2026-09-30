"""Render every delivered PDF page and make readable contact sheets for inspection."""
import math
import fitz
from PIL import Image,ImageOps,ImageDraw
from common import ROOT,dump,sha256

def main():
    dest=ROOT/'paper/qa';dest.mkdir(exist_ok=True);records=[]
    for name in ['manuscript','supplementary','response','cover_letter','manuscript_numbered','manuscript_marked_local_baseline']:
        pdf=ROOT/'paper'/(name+'.pdf');doc=fitz.open(pdf);pages=[];text=[];outside=[]
        for i,page in enumerate(doc):
            pix=page.get_pixmap(matrix=fitz.Matrix(1.6,1.6),alpha=False);path=dest/f'{name}_page_{i+1:02d}.png';pix.save(path);pages.append(path);text.append(page.get_text())
            for block in page.get_text('blocks'):
                if block[0]<5 or block[1]<5 or block[2]>page.rect.width-5 or block[3]>page.rect.height-5:outside.append(dict(page=i+1,bbox=list(block[:4]),text=block[4][:90]))
        contacts=[]
        for k in range(0,len(pages),6):
            subset=pages[k:k+6];thumbs=[]
            for j,path in enumerate(subset):
                a=Image.open(path).convert('RGB');a.thumbnail((490,694));b=Image.new('RGB',(510,726),'#d9dde1');b.paste(a,((510-a.width)//2,24));ImageDraw.Draw(b).text((12,6),f'{name} / page {k+j+1}',fill='black');thumbs.append(b)
            sheet=Image.new('RGB',(1530,726*math.ceil(len(thumbs)/3)),'white')
            for j,a in enumerate(thumbs):sheet.paste(a,((j%3)*510,(j//3)*726))
            q=dest/f'{name}_contact_{k//6+1}.png';sheet.save(q);contacts.append(str(q.relative_to(ROOT)))
        (dest/(name+'_extracted.txt')).write_text('\n\n'.join(text),encoding='utf-8')
        records.append(dict(document=name,pages=len(doc),words_approximate=sum(len(t.split()) for t in text),pdf_sha256=sha256(pdf),page_images=[str(p.relative_to(ROOT)) for p in pages],contact_sheets=contacts,text_near_page_edge=outside));doc.close()
    dump(dest/'render_inventory.json',records)
    for r in records:print(r['document'],r['pages'],'pages;',len(r['text_near_page_edge']),'edge findings',flush=True)

if __name__=='__main__':main()
