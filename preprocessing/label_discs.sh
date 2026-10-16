#!/usr/bin/env bash
# Manual C2-C3 and C7-T1 disc labels for every session that has none. The only interactive step.
#   BIDS_DIR=/path/to/bids ./label_discs.sh
set -euo pipefail

BIDS_DIR="${BIDS_DIR:?set BIDS_DIR}"
DERIV_DIR="${DERIV_DIR:-${BIDS_DIR}/derivatives/sct}"

for d in "$BIDS_DIR"/sub-*/ses-*; do
  sub=$(basename "$(dirname "$d")"); ses=$(basename "$d")
  t2w="$d/anat/${sub}_${ses}_run-1_T2w.nii.gz"
  out="$DERIV_DIR/$sub/$ses/anat"
  [ -f "$t2w" ] && [ ! -f "$out/labels_disc.nii.gz" ] || continue
  mkdir -p "$out"
  [ -f "$out/T2w.nii.gz" ] || cp "$t2w" "$out/T2w.nii.gz"
  echo "$sub $ses"
  ( cd "$out" && sct_label_utils -i T2w.nii.gz -create-viewer 3,8 -o labels_disc.nii.gz \
      -msg "Click the C2-C3 disc, then the C7-T1 disc" )
done
