"""Read-only documentation/runtime-syntax checks; never runs subject processing."""
from pathlib import Path
import argparse
import ast
import hashlib
import json
import re
import subprocess
import sys
from urllib.parse import unquote, urlsplit
import xml.etree.ElementTree as ET

from PIL import Image


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--base-root', type=Path, help='Read-only fallback for resources not copied into a staged overlay')
    parser.add_argument('--tans-root', type=Path, help='Optionally check the real companion driver CLI contract')
    parser.add_argument('--portable-config', action='store_true', help='Require public config to omit personal installation paths')
    args = parser.parse_args()
    root = args.root.resolve()
    base = (args.base_root or root).resolve()
    errors = []
    def require(ok, message):
        if not ok:
            errors.append(message)
    def existing(rel):
        return root/rel if (root/rel).exists() else base/rel
    documents = [root/'README.md', *sorted((root/'docs').rglob('*.md'))]
    links = 0
    for document in documents:
        content = document.read_text()
        hrefs = re.findall(r'!?\[[^\]]*\]\(([^)]+)\)', content)
        hrefs += re.findall(r'(?:href|src)="([^"]+)"', content)
        for href in hrefs:
            parts = urlsplit(href.strip('<>'))
            if parts.scheme or parts.netloc or not parts.path:
                continue
            target = (document.parent/unquote(parts.path)).resolve()
            try:
                rel = target.relative_to(root)
            except ValueError:
                errors.append(f'Link escapes repository: {document.relative_to(root)} -> {href}')
                continue
            path = existing(rel)
            require(path.is_file(), f'Broken link: {document.relative_to(root)} -> {href}')
            if path.is_file() and path.suffix == '.md' and parts.fragment:
                titles = re.findall(r'^#{1,6}\s+(.+)$', path.read_text(), re.M)
                slugs = {re.sub(r'[^\w\- ]', '', t.lower()).replace(' ', '-') for t in titles}
                require(unquote(parts.fragment) in slugs, f'Broken anchor: {href}')
            links += 1
    for relative in ('README.md','config/tans_module.example.sh','config/tans_study.example.yaml',
                     'docs/tans-integration.md','docs/output-gallery.md'):
        require('/home/charleslynch' not in (root/relative).read_text(), f'Personal path in public guide: {relative}')
    if args.portable_config:
        for config in (root/'config').iterdir():
            if config.is_file():
                require(not re.search(r'/home/|/media/|/Users/|charleslynch', config.read_text()),
                        f'Nonportable config: {config.name}')
    for relative in ('modules/mefmri_func_tans.sh','modules/mefmri_func_qc_coreg_sdc.sh',
                     'modules/mefmri_func_qc_motion_grayplots.sh','lib/mefmri_qc_coreg_sdc_report.py',
                     'lib/mefmri_qc_motion_grayplots.py','lib/pfm_ridge_fusion_review.py'):
        require((root/relative).is_file(), f'Missing advertised helper: {relative}')
    shell_files = [p for folder in ('bin','modules','config') for p in (root/folder).glob('*.sh')]
    for path in shell_files:
        result = subprocess.run(['bash','-n',str(path)], capture_output=True, text=True)
        require(result.returncode == 0, f'Shell syntax: {path.name}: {result.stderr}')
    python_files = [p for folder in ('lib','tests','tools') for p in (root/folder).glob('*.py')]
    for path in python_files:
        try:
            ast.parse(path.read_text())
        except SyntaxError as error:
            errors.append(f'Python syntax: {path.name}: {error}')
    manifest = json.loads((root/'docs/assets/manifest.json').read_text())
    names = set()
    for row in manifest['assets']:
        name = row['file']; names.add(name)
        path = root/'docs/assets'/name
        require(path.is_file(), f'Missing asset: {name}')
        if not path.is_file():
            continue
        require(hashlib.sha256(path.read_bytes()).hexdigest() == row['sha256'], f'Changed asset: {name}')
        require(path.stat().st_size < 25*1024**2, f'Asset exceeds 25 MB: {name}')
        if path.suffix in ('.png','.gif'):
            with Image.open(path) as img:
                require(min(img.size)>=350, f'Low-resolution asset: {name}')
                require(not re.search(r'/home/|/tmp/|BioSelect-|CIR002', json.dumps(img.info, default=str)),
                        f'Private path/identifier in image metadata: {name}')
                for frame in range(getattr(img, 'n_frames', 1)):
                    img.seek(frame); img.load()
        elif path.suffix == '.svg':
            ET.parse(path)
        elif path.suffix == '.mp4':
            probe = subprocess.run(['ffprobe','-v','error','-show_entries','stream=width,height','-of','json',str(path)],
                                   capture_output=True, text=True)
            require(probe.returncode == 0, f'Movie decoding metadata: {name}')
    require({p.name for p in (root/'docs/assets').iterdir()} == names | {'manifest.json','README.md'},
            'Unreviewed file in gallery')
    wrapper_options = set(re.findall(r'--[a-z][a-z0-9-]*', (root/'modules/mefmri_func_tans.sh').read_text()))
    if args.tans_root:
        tree = ast.parse((args.tans_root/'pipeline/run_pipeline.py').read_text())
        available = {arg.value for node in ast.walk(tree) if isinstance(node, ast.Call)
                     and isinstance(node.func, ast.Attribute) and node.func.attr == 'add_argument'
                     for arg in node.args if isinstance(arg, ast.Constant) and isinstance(arg.value, str)
                     and arg.value.startswith('--')}
        require(wrapper_options <= available, f'TANS flags absent from companion CLI: {sorted(wrapper_options-available)}')
    print(json.dumps({'status':'failed' if errors else 'passed','documents':len(documents),
        'local_links':links,'shell_files':len(shell_files),'python_files':len(python_files),
        'gallery_assets':len(names),'tans_flags_checked':len(wrapper_options) if args.tans_root else 0,
        'errors':errors,'scope':'Static/read-only; no subject processing or clinical validation.'}, indent=2))
    return bool(errors)

if __name__ == '__main__':
    sys.exit(main())
