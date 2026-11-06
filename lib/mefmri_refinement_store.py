"""Content-addressed reference artifacts; no subject data or pointers are edited."""
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import tempfile


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()


@contextmanager
def publication_lock(root):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    lock = root / '.publish.lock'
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as exc:
        raise RuntimeError('Reference publication is locked; inspect the owning job before removing a stale lock') from exc
    try:
        with os.fdopen(fd, 'w') as stream:
            stream.write(str(os.getpid()) + '\n')
        yield
    finally:
        lock.unlink()


def verify(root, version, compatibility=None):
    if not re.fullmatch('[0-9a-f]{64}', version):
        raise ValueError('An explicit 64-character reference version is required')
    folder = Path(root) / version
    if folder.is_symlink():
        raise ValueError('Reference version directory cannot be a symlink')
    manifest = json.loads((folder / 'manifest.json').read_text())
    if hashlib.sha256(canonical(manifest)).hexdigest() != version:
        raise ValueError('Reference manifest has changed')
    if manifest.get('schema') != 1:
        raise ValueError('Unsupported reference schema')
    for name, checksum in manifest['artifacts'].items():
        if Path(name).name != name or name in ('', '.', '..', 'manifest.json'):
            raise ValueError('Invalid reference artifact name')
        path = folder / name
        if path.is_symlink() or not path.is_file() or digest(path) != checksum:
            raise ValueError('Missing or changed reference artifact: ' + name)
    if compatibility is not None and manifest['compatibility'] != compatibility:
        raise ValueError('Incoming run is incompatible with the frozen reference')
    return manifest


def publish(root, artifacts, members, compatibility, settings):
    """Copy artifacts then atomically publish; identity includes complete provenance.

    members must contain source hashes/transform and motion-target provenance;
    callers construct those from validated run manifests, never path existence.
    """
    if not artifacts or not members:
        raise ValueError('Reference artifacts and explicit membership are required')
    if len({m['id'] for m in members}) != len(members):
        raise ValueError('Duplicate reference members')
    for name in artifacts:
        if Path(name).name != name or name in ('', '.', '..', 'manifest.json'):
            raise ValueError('Invalid artifact filename')
    root = Path(root)
    with publication_lock(root):
        stage = Path(tempfile.mkdtemp(prefix='.building-', dir=root))
        try:
            for name, source in artifacts.items():
                shutil.copyfile(source, stage / name)
            manifest = dict(schema=1, members=sorted(members, key=lambda m: m['id']),
                            compatibility=compatibility, settings=settings,
                            artifacts={name: digest(stage / name) for name in sorted(artifacts)})
            version = hashlib.sha256(canonical(manifest)).hexdigest()
            (stage / 'manifest.json').write_bytes(canonical(manifest) + b'\n')
            destination = root / version
            if destination.exists():
                verify(root, version, compatibility)
            else:
                stage.rename(destination)
            return version
        finally:
            if stage.exists():
                shutil.rmtree(stage)


def require_run_binding(saved, expected):
    """Reject reuse of downstream BOLD generated under a different alignment."""
    if canonical(saved) != canonical(expected):
        raise ValueError('Alignment inputs/version changed: regenerate spatial BOLD and its downstream products')
