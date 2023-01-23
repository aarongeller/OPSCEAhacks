# OPSCEAhacks

## Overview
My mods to the [UCSF OPSCEA project](https://github.com/Kleen-Lab/OPSCEA) (March 2022 version).  Main functions provided:
1. Using [Brainstorm](https://neuroimage.usc.edu/brainstorm/Introduction) for electrode localization;
2. A function, `OPSCEA_struct`, that exports slice images (without EEG heatmaps, including both coronal and axial/sagittal view) and compiles a pdf document of the reconstruction.

## Workflow
1. Run freesurfer.
2. In Brainstorm,
 * Load MRI and co-register CT
 * Mark electrodes
3. Run `create_opscea_subj.py`
4. Edit OPSCEAparams.xls
5. Run `OPSCEA.m` or `OPSCEA_struct.m`.