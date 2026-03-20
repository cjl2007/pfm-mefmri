#!/bin/bash
# CJL; (cjl2007@med.cornell.edu)

Subject=$1
StudyFolder=$2
Subdir="$StudyFolder"/"$Subject"
MEDIR=$3
StartSession=$4
FuncDirName="${FUNC_DIRNAME:-rest}"
FuncFilePrefix="${FUNC_FILE_PREFIX:-Rest}"

# python runtime override (default: python3 on PATH)
: "${MGTR_PYTHON:=python3}"
MGTR_PY_SCRIPT="$MEDIR/lib/mgtr_volume.py"
if [ ! -f "$MGTR_PY_SCRIPT" ]; then
	echo "ERROR: missing MGTR Python script: $MGTR_PY_SCRIPT"
	exit 2
fi

# count the number of sessions
sessions=("$Subdir"/func/"$FuncDirName"/session_*)
sessions=$(seq $StartSession 1 "${#sessions[@]}")

# sweep the sessions;
for s in $sessions ; do

	# count number of runs for this session;
	runs=("$Subdir"/func/"$FuncDirName"/session_"$s"/run_*)
	runs=$(seq 1 1 "${#runs[@]}" )

	# Iterate over runs.
	for r in $runs ; do

		Input="$Subdir/func/$FuncDirName/session_$s/run_$r/${FuncFilePrefix}_OCME+MEICA.nii.gz"
		Output_MGTR="$Subdir/func/$FuncDirName/session_$s/run_$r/${FuncFilePrefix}_OCME+MEICA+MGTR"
		Output_Betas="$Subdir/func/$FuncDirName/session_$s/run_$r/${FuncFilePrefix}_OCME+MEICA+MGTR_Betas"

		if [ ! -f "$Input" ]; then
			echo "ERROR: missing MGTR input: $Input"
			exit 2
		fi

		"$MGTR_PYTHON" "$MGTR_PY_SCRIPT" \
			--subdir "$Subdir" \
			--input "$Input" \
			--output-mgtr-base "$Output_MGTR" \
			--output-betas-base "$Output_Betas"

	done
	
done
