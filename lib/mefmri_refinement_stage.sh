# Wrapper integration. Baseline pointers remain intact.
refinement_stage() {
  local action="$1" resume="${2:-0}" result mode="$DISTORTION_CORRECTION_MODE"
  case "$mode" in fieldmap|gre_fieldmap) mode=direct_b0 ;; esac
  local args=("$action" --subject-dir "$SubjectDir" --func-dir "$FUNC_DIRNAME"
    --xfms-dir "$FUNC_XFMS_DIRNAME" --mode "$mode" --space "$AtlasSpace" --grid "$AtlasTemplate"
    --policy "$COREG_REFERENCE_POLICY" --strategy "$COREG_ALIGNMENT_STRATEGY" --start-session "$START_SESSION")
  [[ -z "${COREG_ALIGNMENT_BUNDLE:-}" ]] || args+=(--bundle "$COREG_ALIGNMENT_BUNDLE")
  [[ -z "${COREG_REFERENCE_VERSION:-}" ]] || args+=(--version "$COREG_REFERENCE_VERSION")
  [[ -z "${COREG_REFERENCE_MEMBERS:-}" ]] || args+=(--members "$COREG_REFERENCE_MEMBERS")
  [[ -z "${COREG_REFERENCE_ROOT:-}" ]] || args+=(--reference-root "$COREG_REFERENCE_ROOT")
  [[ "$resume" != 1 ]] || args+=(--resume)
  result=$("${PIPELINE_PYTHON:-python3}" "$MEDIR/lib/mefmri_refinement_stage.py" "${args[@]}") || return 1
  if [[ "$action" == prepare && "$COREG_ALIGNMENT_STRATEGY" == subject_refined ]]; then
    local fields=()
    mapfile -t fields <<< "$result"
    [[ ${#fields[@]} -eq 3 ]] || { echo '[ERROR] Invalid alignment-stage response' >&2; return 1; }
    export COREG_ALIGNMENT_BUNDLE="${fields[0]}" COREG_REFERENCE_VERSION="${fields[1]}"
    export FUNC_XFMS_DIRNAME="${fields[2]}"
    echo "[coreg] Subject-refined bundle: $COREG_ALIGNMENT_BUNDLE"
    echo "[coreg] Reference version: $COREG_REFERENCE_VERSION"
    if [[ -n "${RUN_META_FILE:-}" ]]; then
      {
        echo "coreg_reference_version_resolved=$COREG_REFERENCE_VERSION"
        echo "coreg_alignment_bundle_resolved=$COREG_ALIGNMENT_BUNDLE"
        echo "coreg_xfms_namespace_resolved=$FUNC_XFMS_DIRNAME"
      } >> "$RUN_META_FILE"
    fi
  elif [[ "$action" == reserve ]]; then
    local fields=()
    mapfile -t fields <<< "$result"
    [[ ${#fields[@]} -eq 2 ]] || return 1
    export FUNC_XFMS_DIRNAME="${fields[0]}" COREG_REFERENCE_ROOT="${fields[1]}"
    echo "[coreg] Separate baseline transform namespace: $FUNC_XFMS_DIRNAME"
    if [[ -n "${RUN_META_FILE:-}" ]]; then
      {
        echo "coreg_reference_root_resolved=$COREG_REFERENCE_ROOT"
        echo "coreg_xfms_namespace_reserved=$FUNC_XFMS_DIRNAME"
      } >> "$RUN_META_FILE"
    fi
  fi
}
