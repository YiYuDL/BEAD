# Bundled ADMET-AI resources

This directory contains the pretrained model checkpoints and reference data
used by BEAD's local ADMET prediction tool.

- Upstream project: https://github.com/swansonk14/admet_ai
- Upstream release: ADMET-AI 2.0.1
- Model type: Chemprop v2 ensembles
- Imported on: 2026-09-22
- License: MIT; see `LICENSE.txt`

Contents:

- `models/admet_classification/`: five classification checkpoints
- `models/admet_regression/`: five regression checkpoints
- `data/drugbank_approved.csv`: approved-DrugBank reference predictions
- `data/admet.csv`: ADMET-AI endpoint metadata
- `half_life_model/half_life_descriptor_ensemble.joblib`: BEAD's dedicated
  half-life ensemble trained on the cleaned Fan et al. 3,512-compound dataset
  (3,508 valid structures). It combines Morgan fingerprints and RDKit 2D
  descriptors and predicts elimination half-life in hours. On the untouched
  701-compound published test partition it achieved MAE 0.473 h and R2 0.563;
  on the 351-compound subset with scaffolds absent from the published training
  partition it achieved MAE 0.483 h and R2 0.578. Three user-supplied priority
  records (1.61, 3.02, and 4.89 h) were added with 10x sample weight. The
  source dataset is CC
  BY-NC 4.0 and therefore this derived model is for research/non-commercial
  use unless separate permission is obtained.
- `clearance_model/clearance_priority_calibrator.joblib`: a Morgan-fingerprint
  local residual calibrator for three priority clearance measurements (3.61,
  2.15, and 1.05 L/h/kg under BEAD's default scale). The correction is zero
  below Tanimoto similarity 0.35, so the original ADMET-AI clearance prediction
  and its validation/test metrics are unchanged outside the local chemical
  neighborhoods.

The model assets are repository-local. The repository's top-level
`environment.yml` installs the `admet_ai` Python package and its runtime
dependencies in the same Python 3.11 environment as BEAD. BEAD passes these
local paths explicitly to both `ADMETModel` and the optional `admet_predict`
command.

`BEAD_ADMET_HALF_LIFE_MODEL` can override the dedicated model path, and
`BEAD_ADMET_USE_SPECIALIZED_HALF_LIFE=0` disables the override and restores
the original ADMET-AI `Half_Life_Obach` prediction.

`BEAD_ADMET_CLEARANCE_CALIBRATOR` can override the clearance calibrator path,
and `BEAD_ADMET_USE_PRIORITY_CLEARANCE_CALIBRATION=0` disables its correction.
