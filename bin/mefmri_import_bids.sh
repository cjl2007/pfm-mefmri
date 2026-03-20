#!/usr/bin/env bash
# Convert BIDS inputs into RevisedMe-fMRIPipeline expected raw folder layout.
#
# Usage:
#   mefmri_import_bids.sh <BIDS_ROOT> <SUBJECT> <OUT_SUBJECT_DIR> [options]
# Example:
#   mefmri_import_bids.sh /data/bids 06 /data/study/ME06 --task rest --mode symlink --overwrite

set -euo pipefail

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  MEDIR="$(cd "$SCRIPT_DIR/.." && pwd)"
  python3 "$MEDIR/lib/mefmri_bids_import.py" --help
  exit 0
fi

if [[ "$#" -lt 3 ]]; then
  echo "Usage: $0 <BIDS_ROOT> <SUBJECT> <OUT_SUBJECT_DIR> [--task rest] [--func-dirname rest] [--func-prefix Rest] [--mode symlink|copy] [--overwrite]" >&2
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MEDIR="$(cd "$SCRIPT_DIR/.." && pwd)"

python3 "$MEDIR/lib/mefmri_bids_import.py" "$@"
