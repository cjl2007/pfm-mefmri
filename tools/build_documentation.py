#!/usr/bin/env python3
"""Build portable HTML guides using the shared documentation stylesheet."""
from pathlib import Path
import argparse,html,json,os,re,subprocess
from urllib.parse import unquote,urlsplit

def build(root,output,overlay=None):
    root=Path(root).resolve();output=Path(output).resolve();overlay=Path(overlay).resolve() if overlay else None
    sources={p.relative_to(root):p for p in (root/'docs').rglob('*.md') if not any('backup' in s.lower() for s in p.parts)}
    sources[Path('README.md')]=root/'README.md'
    if overlay:
        if (overlay/'README.md').is_file():sources[Path('README.md')]=overlay/'README.md'
        for p in (overlay/'docs').rglob('*.md'):sources[p.relative_to(overlay)]=p
    pages={rel:Path('docs/index.html') if rel==Path('README.md') else rel.with_suffix('.html') for rel in sources}
    links=[]
    for rel,source in sources.items():
        dest=output/pages[rel];dest.parent.mkdir(parents=True,exist_ok=True)
        body=subprocess.run(['pandoc','--from=gfm','--to=html5','--wrap=none',str(source)],check=True,text=True,capture_output=True).stdout
        def url(match):
            attr,value=match.group(1),html.unescape(match.group(2));u=urlsplit(value)
            if not u.path:return match.group(0)
            if u.scheme:
                prefix='https://raw.githubusercontent.com/cjl2007/pfm-mefmri/main/'
                if not value.startswith(prefix):return match.group(0)
                target=Path(unquote(u.path.split('/main/',1)[1]))
            else:target=Path(os.path.normpath(str(rel.parent/unquote(u.path))))
            actual=pages.get(target,target)
            exists=target in pages or (overlay and (overlay/target).exists()) or (root/target).exists()
            links.append(dict(page=str(rel),target=str(target),exists=bool(exists)))
            location=os.path.relpath(actual,pages[rel].parent)
            if u.fragment:location+='#'+u.fragment
            return attr+'="'+html.escape(location,quote=True)+'"'
        def image_box(match):
            tag=match.group(0);src=html.unescape(re.search(r'src="([^"]+)"',tag).group(1))
            u=urlsplit(src)
            if u.scheme:return tag
            target=Path(os.path.normpath(str(rel.parent/unquote(u.path))))
            source=(overlay/target) if overlay and (overlay/target).is_file() else root/target
            if source.suffix.lower() not in ['.png','.gif','.jpg','.jpeg']:return tag
            from PIL import Image
            with Image.open(source) as im:w,h=im.size
            return tag[:-1]+' width="'+str(w)+'" height="'+str(h)+'" style="height:auto;aspect-ratio:'+str(w)+' / '+str(h)+'"'+tag[-1:]
        body=re.sub(r'<img\b[^>]*src="[^"]+"[^>]*>',image_box,body)
        body=re.sub(r'(src|href|poster)="([^"]+)"',url,body)
        css=os.path.relpath(Path('docs/assets/docs.css'),pages[rel].parent)
        home=os.path.relpath(Path('docs/index.html'),pages[rel].parent)
        guides=os.path.relpath(Path('docs/README.html'),pages[rel].parent)
        heading=re.search(r'<h1[^>]*>(.*?)</h1>',body)
        title=re.sub('<[^>]+>','',heading.group(1)) if heading else 'pfm-mefmri'
        dest.write_text('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+title+'</title><link rel="stylesheet" href="'+css+'"></head><body><main><nav><a href="'+home+'">pfm-mefmri</a><a href="'+guides+'">Guides</a></nav>'+body+'</main></body></html>')
    # Retain the existing short link to the TANS guide.
    (output/'docs/tans.html').write_text('<!doctype html><meta charset="utf-8"><meta http-equiv="refresh" content="0;url=tans/README.html"><a href="tans/README.html">Individualized TMS</a>')
    report=dict(pages=len(pages),links=len(links),missing=[r for r in links if not r['exists']])
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);p.add_argument('--output-root',type=Path);p.add_argument('--overlay',type=Path);a=p.parse_args()
    result=build(a.root,a.output_root or a.root,a.overlay);print(json.dumps(result))
    if result['missing']:raise SystemExit(1)
