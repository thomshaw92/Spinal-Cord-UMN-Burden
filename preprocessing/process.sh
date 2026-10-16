#!/usr/bin/env bash
# Spinal Cord Toolbox 6.5 processing of one session, or of every session in the BIDS folder.
#   BIDS_DIR=/path/to/bids ./process.sh [sub-XX ses-XX]
# Run label_discs.sh first: registration uses the manual C2-C3 and C7-T1 disc labels.
set -uo pipefail

BIDS_DIR="${BIDS_DIR:?set BIDS_DIR}"
DERIV_DIR="${DERIV_DIR:-${BIDS_DIR}/derivatives/sct}"
PAM50="${SCT_DIR:-$HOME/spinalcordtoolbox}/data/PAM50/template"
LEVELS="1:9"

t2w() {
  cd "$1/anat"
  [ -f T2w.nii.gz ] || cp "$2/anat/$3_run-1_T2w.nii.gz" T2w.nii.gz
  [ -f T2w_seg.nii.gz ] || sct_deepseg_sc -i T2w.nii.gz -c t2 -o T2w_seg.nii.gz -qc qc
  [ -f labels_disc.nii.gz ] || { echo "no disc labels for $3, run label_discs.sh"; return 1; }
  sct_register_to_template -i T2w.nii.gz -s T2w_seg.nii.gz -ldisc labels_disc.nii.gz -c t2 -qc qc
  sct_warp_template -d T2w.nii.gz -w warp_template2anat.nii.gz -a 1 -qc qc
  sct_process_segmentation -i T2w_seg.nii.gz -vert "$LEVELS" -perlevel 1 \
    -vertfile label/template/PAM50_levels.nii.gz -o csa_T2w_perlevel.csv -qc qc
}

register() {
  sct_register_multimodal -i "$PAM50/$1" -iseg "$PAM50/PAM50_cord.nii.gz" -d "$2" -dseg "$3" \
    -initwarp "$4/warp_template2anat.nii.gz" -initwarpinv "$4/warp_anat2template.nii.gz" \
    -owarp "warp_template2$5.nii.gz" -owarpinv "warp_${5}2template.nii.gz" -param "$6" -qc qc
  rm -rf "label_$5"
  sct_warp_template -d "$2" -w "warp_template2$5.nii.gz" -a 1 -ofolder "label_$5" -qc qc
}

t2star() {
  cd "$1/anat"
  [ -f T2starw.nii.gz ] || cp "$2/anat/$3_run-1_T2starw.nii.gz" T2starw.nii.gz || return 1
  [ -f T2starw_seg.nii.gz ] || sct_deepseg_sc -i T2starw.nii.gz -c t2s -o T2starw_seg.nii.gz -qc qc
  [ -f T2starw_gmseg.nii.gz ] || sct_deepseg_gm -i T2starw.nii.gz -o T2starw_gmseg.nii.gz -qc qc
  [ -f T2starw_wmseg.nii.gz ] || sct_maths -i T2starw_seg.nii.gz -sub T2starw_gmseg.nii.gz -o T2starw_wmseg.nii.gz
  register PAM50_t2s.nii.gz T2starw.nii.gz T2starw_seg.nii.gz . T2s \
    "step=1,type=seg,algo=centermass:step=2,type=seg,algo=bsplinesyn,slicewise=1,iter=3:step=3,type=im,algo=syn,slicewise=1,iter=1,metric=CC"
  for s in gma:T2starw_gmseg cord_csa:T2starw_seg wm_csa:T2starw_wmseg; do
    sct_process_segmentation -i "${s#*:}.nii.gz" -vertfile label_T2s/template/PAM50_levels.nii.gz \
      -vert "$LEVELS" -perlevel 1 -o "${s%%:*}_T2starw_perlevel.csv" -qc qc
  done
}

dwi() {
  cd "$1/dwi"
  local raw="$2/dwi/$3_run-1_dwi"
  [ -f dwi.nii.gz ] || cp "$raw.nii.gz" dwi.nii.gz || return 1
  [ -f bvals.txt ] || cp "$raw.bval" bvals.txt
  [ -f bvecs.txt ] || cp "$raw.bvec" bvecs.txt
  [ -f dwi_dwi_mean.nii.gz ] || sct_dmri_separate_b0_and_dwi -i dwi.nii.gz -bvec bvecs.txt -a 1
  [ -f dwi_seg_premoco.nii.gz ] || sct_deepseg_sc -i dwi_dwi_mean.nii.gz -c dwi -o dwi_seg_premoco.nii.gz -qc qc
  [ -f dwi_mask.nii.gz ] || sct_create_mask -i dwi_dwi_mean.nii.gz -p centerline,dwi_seg_premoco.nii.gz -size 35mm -o dwi_mask.nii.gz
  if [ ! -f dwi_moco.nii.gz ]; then
    # MeanSquares when the default mutual-information metric fails on low-contrast volumes
    sct_dmri_moco -i dwi.nii.gz -bvec bvecs.txt -bval bvals.txt -m dwi_mask.nii.gz -qc qc -qc-seg dwi_seg_premoco.nii.gz ||
      sct_dmri_moco -i dwi.nii.gz -bvec bvecs.txt -bval bvals.txt -m dwi_mask.nii.gz -param metric=MeanSquares -qc qc -qc-seg dwi_seg_premoco.nii.gz
  fi
  [ -f dwi_seg.nii.gz ] || sct_deepseg_sc -i dwi_moco_dwi_mean.nii.gz -c dwi -o dwi_seg.nii.gz -qc qc
  [ -f dti_FA.nii.gz ] || sct_dmri_compute_dti -i dwi_moco.nii.gz -bval bvals.txt -bvec bvecs.txt -m dwi_seg.nii.gz -o dti_ -method standard
  register PAM50_t1.nii.gz dwi_moco_dwi_mean.nii.gz dwi_seg.nii.gz ../anat dwi \
    "step=1,type=seg,algo=centermass:step=2,type=seg,algo=bsplinesyn,slicewise=1,iter=3"
  for m in dti_FA dti_MD dti_AD dti_RD; do
    sct_extract_metric -i "$m.nii.gz" -f label_dwi/atlas -method map -vert "$LEVELS" -perlevel 1 \
      -vertfile label_dwi/template/PAM50_levels.nii.gz -l 4,5,22,23,0,1,2,3,12,13 -o "${m}_tracts_perlevel.csv"
  done
}

mt() {
  cd "$1/anat"
  [ -f MTon.nii.gz ] || cp "$2/anat/$3_acq-MTon_run-1_MTS.nii.gz" MTon.nii.gz || return 1
  [ -f MToff.nii.gz ] || cp "$2/anat/$3_acq-MToff_run-1_MTS.nii.gz" MToff.nii.gz || return 1
  [ -f MTon_seg.nii.gz ] || sct_deepseg_sc -i MTon.nii.gz -c t2 -o MTon_seg.nii.gz -qc qc
  register PAM50_t2.nii.gz MTon.nii.gz MTon_seg.nii.gz . MT \
    "step=1,type=seg,algo=centermass:step=2,type=seg,algo=bsplinesyn,slicewise=1,iter=3"
  [ -f MToff_reg.nii.gz ] || sct_register_multimodal -i MToff.nii.gz -d MTon.nii.gz -dseg MTon_seg.nii.gz \
    -param step=1,type=im,algo=rigid,metric=CC
  [ -f mtr.nii.gz ] || sct_compute_mtr -mt0 MToff_reg.nii.gz -mt1 MTon.nii.gz -o mtr.nii.gz
  sct_extract_metric -i mtr.nii.gz -f label_MT/atlas -method map -vert "$LEVELS" -perlevel 1 \
    -vertfile label_MT/template/PAM50_levels.nii.gz -l 4,5,22,23,0,1,2,3,12,13,30,31,34,35 -o mtr_tracts_perlevel.csv
}

session() {
  local sub="$1" ses="$2" src="$BIDS_DIR/$1/$2" out="$DERIV_DIR/$1/$2"
  mkdir -p "$out/anat" "$out/dwi"
  ( t2w "$out" "$src" "${sub}_${ses}" ) || { echo "T2w failed: $sub $ses"; return; }
  for stage in t2star dwi mt; do
    ( "$stage" "$out" "$src" "${sub}_${ses}" ) || echo "$stage failed or missing: $sub $ses"
  done
}

if [ $# -eq 2 ]; then
  session "$1" "$2"
else
  for d in "$BIDS_DIR"/sub-*/ses-*; do
    session "$(basename "$(dirname "$d")")" "$(basename "$d")"
  done
fi
