#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  build_container_image.sh [--tag IMAGE_TAG] [--base BASE_IMAGE]

Defaults:
  IMAGE_TAG=mefmri-pipeline:runtime
  BASE_IMAGE=nipreps/fmriprep:24.1.1
USAGE
}

IMAGE_TAG="mefmri-pipeline:runtime"
BASE_IMAGE="nipreps/fmriprep:24.1.1"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --tag)
      IMAGE_TAG="${2:?missing value for --tag}"
      shift 2
      ;;
    --base)
      BASE_IMAGE="${2:?missing value for --base}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: unknown argument: $1" >&2
      usage
      exit 2
      ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

cd "$REPO_ROOT"
docker build \
  -f docker/Dockerfile.runtime \
  --build-arg BASE_IMAGE="$BASE_IMAGE" \
  -t "$IMAGE_TAG" \
  .

echo "Built image: $IMAGE_TAG"
