"""Bind downstream completion to the precise spatial BOLD generation."""
import json
import os
from pathlib import Path
import tempfile
from mefmri_refinement_store import digest


STAGES = {'denoise', 'mgtr', 'surface', 'concat'}


class AlignmentRecordNeedsRegeneration(ValueError):
    """Expected status when a stage has not run for the current BOLD."""


class MissingAlignmentRecord(AlignmentRecordNeedsRegeneration):
    """No completion record exists yet for this stage."""


def record_path(run, stage):
    if stage not in STAGES:
        raise ValueError('Unknown downstream alignment stage')
    return Path(run) / (stage.capitalize() + 'Alignment.json')


def check(run, stage):
    run = Path(run)
    bold = run / 'BOLDAlignment.json'
    if not bold.exists():
        return  # Legacy independent runs have no alignment lifecycle records.
    path = record_path(run, stage)
    if not path.is_file():
        raise MissingAlignmentRecord('Missing ' + stage + ' alignment record; regenerate that stage from current BOLD: ' + str(run))
    saved = json.loads(path.read_text())
    if saved.get('bold_record_sha256') != digest(bold):
        raise AlignmentRecordNeedsRegeneration('Stale ' + stage + ' outputs after BOLD/alignment change: ' + str(run))


def stamp(run, stage):
    run = Path(run)
    bold = run / 'BOLDAlignment.json'
    if not bold.exists():
        return
    record = dict(schema=1, stage=stage, bold_record_sha256=digest(bold))
    fd, name = tempfile.mkstemp(prefix='.downstream-alignment-', dir=run)
    try:
        with os.fdopen(fd, 'w') as stream:
            json.dump(record, stream, sort_keys=True)
            stream.write('\n')
        os.replace(name, record_path(run, stage))
    finally:
        if Path(name).exists():
            Path(name).unlink()


def process(subject, func_dir, start_session, stage, action):
    if action not in ('check', 'stamp'):
        raise ValueError('Invalid downstream action')
    for run in (Path(subject) / 'func' / func_dir).glob('session_*/run_*'):
        try:
            session = int(run.parent.name.split('session_', 1)[-1])
        except ValueError:
            continue
        if session >= start_session:
            (check if action == 'check' else stamp)(run, stage)


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['check', 'stamp'])
    p.add_argument('--subject-dir', required=True)
    p.add_argument('--func-dir', required=True)
    p.add_argument('--start-session', type=int, default=1)
    p.add_argument('--stage', choices=sorted(STAGES), required=True)
    p.add_argument('--quiet-missing', action='store_true',
                   help='Suppress missing-record diagnostics when the stage will run next')
    a = p.parse_args()
    try:
        process(a.subject_dir, a.func_dir, a.start_session, a.stage, a.action)
    except AlignmentRecordNeedsRegeneration as exc:
        # Nonzero status is intentional: the denoise wrapper regenerates
        # outputs, while late-stage resumes still stop on missing/stale data.
        # Unexpected exceptions retain their full diagnostic traceback.
        message = None if a.quiet_missing and isinstance(exc, MissingAlignmentRecord) else '[alignment] ' + str(exc) + '\n'
        p.exit(1, message)
