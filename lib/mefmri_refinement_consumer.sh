# Source inside a headmotion run worker. Selection returns data, never shell code.
refinement_select() {
  local pipeline="$1" run_dir="$2" values
  case "${COREG_ALIGNMENT_STRATEGY:-independent}" in
    independent) return 0 ;;
    subject_refined) ;;
    *) echo "[ERROR] Unknown coregistration alignment strategy" >&2; return 1 ;;
  esac
  values=$("${PIPELINE_PYTHON:-python3}" "$pipeline/lib/mefmri_refinement_consumer.py" select \
    --run-dir "$run_dir" --strategy subject_refined \
    --bundle "${COREG_ALIGNMENT_BUNDLE:?Missing alignment bundle}" \
    --version "${COREG_REFERENCE_VERSION:?Missing reference version}") || return 1
  local fields=()
  mapfile -t fields <<< "$values"
  [[ ${#fields[@]} -ge 2 ]] || return 1
  IntermediateCoregTarget="${fields[0]}"
  Intermediate2ACPCWarp="${fields[1]}"
  Intermediate2ACPCMat="${fields[2]:-}"
  if [[ "${FUNC_NOFIELDMAP_MODE:-0}" == "1" && ! -f "$Intermediate2ACPCMat" ]]; then
    echo "[ERROR] Refined no-fieldmap affine is missing" >&2
    return 1
  fi
}

refinement_mark() {
  local pipeline="$1" run_dir="$2" prefix="$3"
  local args=(--run-dir "$run_dir" --strategy "${COREG_ALIGNMENT_STRATEGY:-independent}" --prefix "$prefix")
  if [[ "${COREG_ALIGNMENT_STRATEGY:-independent}" == "subject_refined" ]]; then
    args+=(--bundle "$COREG_ALIGNMENT_BUNDLE" --version "$COREG_REFERENCE_VERSION")
  fi
  "${PIPELINE_PYTHON:-python3}" "$pipeline/lib/mefmri_refinement_consumer.py" mark "${args[@]}"
}
