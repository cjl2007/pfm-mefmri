"""Reserve separate baseline transform directories before running coregistration."""
from pathlib import Path
import re
import uuid
from mefmri_refinement_store import verify


def reserve(subject_dir, group, policy, version=None, start_session=1, root=None):
    subject_dir = Path(subject_dir).resolve()
    if not re.fullmatch('[A-Za-z0-9_.-]+', group) or group in ('.', '..'):
        raise ValueError('Invalid baseline transform group')
    root = Path(root).resolve() if root else subject_dir / 'func/xfms' / group / 'SubjectRefinement'
    if policy == 'frozen':
        if not version:
            raise ValueError('Frozen processing requires an explicit reference version')
        reference = verify(root / 'references', version)
        if reference['compatibility']['subject'] != subject_dir.name:
            raise ValueError('Frozen reference belongs to a different subject')
        # Coreg's current selection is session-based. Never let an incoming
        # invocation overwrite baseline inputs of calibration runs.
        if any(m['session'] >= start_session for m in reference['members']):
            raise ValueError('Frozen incoming processing must start after calibration sessions; use a saved bundle to resume calibration runs')
    elif policy == 'rebuild':
        if start_session != 1 or version:
            raise ValueError('Rebuild requires START_SESSION=1 and no existing reference version')
    else:
        raise ValueError('Invalid reference policy')
    name = group + '_refinement_' + uuid.uuid4().hex[:12]
    folder = subject_dir / 'func/xfms' / name
    folder.mkdir(parents=True, exist_ok=False)
    root.mkdir(parents=True, exist_ok=True)
    return name, str(root)
