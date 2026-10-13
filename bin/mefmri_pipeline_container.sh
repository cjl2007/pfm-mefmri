#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  mefmri_pipeline_container.sh [--engine docker|apptainer] [--image IMAGE] [--fs-license PATH] <SubjectDir> [ConfigFile]

Environment defaults:
  CONTAINER_ENGINE          docker
  MEFMRI_CONTAINER_IMAGE    mefmri-pipeline:runtime
  FS_LICENSE / FS_LICENSE_FILE (used if --fs-license is not provided)

Examples:
  bash bin/mefmri_pipeline_container.sh /data/study/ME06 config/mefmri_wrapper_config_container.sh
  bash bin/mefmri_pipeline_container.sh --engine apptainer --image mefmri-pipeline.sif --fs-license /path/license.txt /data/study/ME06
USAGE
}

ENGINE="${CONTAINER_ENGINE:-docker}"
IMAGE="${MEFMRI_CONTAINER_IMAGE:-mefmri-pipeline:runtime}"
FS_LICENSE_ARG=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --engine)
      ENGINE="${2:?Missing value for --engine}"
      shift 2
      ;;
    --image)
      IMAGE="${2:?Missing value for --image}"
      shift 2
      ;;
    --fs-license)
      FS_LICENSE_ARG="${2:?Missing value for --fs-license}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      break
      ;;
    -*)
      echo "ERROR: unknown option: $1" >&2
      usage
      exit 2
      ;;
    *)
      break
      ;;
  esac
done

if [[ $# -lt 1 || $# -gt 2 ]]; then
  usage
  exit 2
fi

SubjectDir="$1"
ConfigFileArg="${2:-}"

if [[ "${SubjectDir: -1}" == "/" ]]; then
  SubjectDir="${SubjectDir%?}"
fi
SubjectDir="$(cd "$SubjectDir" && pwd)"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
DEFAULT_CONFIG="$REPO_ROOT/config/mefmri_wrapper_config_container.sh"
ConfigFile="${ConfigFileArg:-${CONFIG_FILE:-$DEFAULT_CONFIG}}"
ConfigFile="$(cd "$(dirname "$ConfigFile")" && pwd)/$(basename "$ConfigFile")"

if [[ ! -d "$SubjectDir" ]]; then
  echo "ERROR: subject directory not found: $SubjectDir" >&2
  exit 2
fi
if [[ ! -f "$ConfigFile" ]]; then
  echo "ERROR: config file not found: $ConfigFile" >&2
  exit 2
fi

FS_LICENSE_HOST="${FS_LICENSE_ARG:-${FS_LICENSE:-${FS_LICENSE_FILE:-}}}"
if [[ -n "$FS_LICENSE_HOST" ]]; then
  FS_LICENSE_HOST="$(cd "$(dirname "$FS_LICENSE_HOST")" && pwd)/$(basename "$FS_LICENSE_HOST")"
  if [[ ! -f "$FS_LICENSE_HOST" ]]; then
    echo "ERROR: FreeSurfer license file not found: $FS_LICENSE_HOST" >&2
    exit 2
  fi
fi

StudyFolder="$(dirname "$SubjectDir")"
SubjectName="$(basename "$SubjectDir")"
ConfigBase="$(basename "$ConfigFile")"
ConfigInContainer=""
ExtraConfigMountDir=""
if [[ "$ConfigFile" == "$REPO_ROOT/"* ]]; then
  ConfigRel="${ConfigFile#"$REPO_ROOT"/}"
  ConfigInContainer="/opt/mefmri/${ConfigRel}"
else
  ExtraConfigMountDir="$(dirname "$ConfigFile")"
  ConfigInContainer="/data/config/$ConfigBase"
fi
FS_LICENSE_CONT=""
ExtraLicenseMountDir=""
if [[ -n "$FS_LICENSE_HOST" ]]; then
  if [[ "$FS_LICENSE_HOST" == "$REPO_ROOT/"* ]]; then
    LicenseRel="${FS_LICENSE_HOST#"$REPO_ROOT"/}"
    FS_LICENSE_CONT="/opt/mefmri/${LicenseRel}"
  else
    ExtraLicenseMountDir="$(dirname "$FS_LICENSE_HOST")"
    FS_LICENSE_CONT="/licenses/$(basename "$FS_LICENSE_HOST")"
  fi
fi

run_docker() {
  local -a cmd=(
    docker run --rm -it
    -v "$REPO_ROOT:/opt/mefmri:ro"
    -v "$StudyFolder:/data/study"
    -e "MEDIR=/opt/mefmri"
  )
  if [[ -n "$ExtraConfigMountDir" ]]; then
    cmd+=(-v "$ExtraConfigMountDir:/data/config:ro")
  fi
  if [[ -n "$FS_LICENSE_HOST" ]]; then
    if [[ -n "$ExtraLicenseMountDir" ]]; then
      cmd+=(-v "$ExtraLicenseMountDir:/licenses:ro")
    fi
    cmd+=(-e "FS_LICENSE=$FS_LICENSE_CONT")
    cmd+=(-e "FS_LICENSE_FILE=$FS_LICENSE_CONT")
  fi
  cmd+=(
    "$IMAGE"
    "/data/study/$SubjectName"
    "$ConfigInContainer"
  )
  "${cmd[@]}"
}

run_apptainer() {
  local -a env_args=(
    "MEDIR=/opt/mefmri"
  )
  if [[ -n "$FS_LICENSE_CONT" ]]; then
    env_args+=("FS_LICENSE=${FS_LICENSE_CONT}" "FS_LICENSE_FILE=${FS_LICENSE_CONT}")
  fi
  local -a cmd=(
    apptainer exec --cleanenv
    -B "$REPO_ROOT:/opt/mefmri"
    -B "$StudyFolder:/data/study"
  )
  if [[ -n "$ExtraConfigMountDir" ]]; then
    cmd+=(-B "$ExtraConfigMountDir:/data/config")
  fi
  if [[ -n "$FS_LICENSE_HOST" ]]; then
    if [[ -n "$ExtraLicenseMountDir" ]]; then
      cmd+=(-B "$ExtraLicenseMountDir:/licenses")
    fi
  fi
  cmd+=(
    "$IMAGE"
    env
    "${env_args[@]}"
    bash /opt/mefmri/bin/mefmri_pipeline.sh
    "/data/study/$SubjectName"
    "$ConfigInContainer"
  )
  "${cmd[@]}"
}

case "$ENGINE" in
  docker)
    command -v docker >/dev/null 2>&1 || { echo "ERROR: docker not found in PATH" >&2; exit 2; }
    run_docker
    ;;
  apptainer|singularity)
    command -v apptainer >/dev/null 2>&1 || { echo "ERROR: apptainer not found in PATH" >&2; exit 2; }
    run_apptainer
    ;;
  *)
    echo "ERROR: unsupported container engine '$ENGINE' (use docker or apptainer)" >&2
    exit 2
    ;;
esac
