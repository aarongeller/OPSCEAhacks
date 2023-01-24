# OPSCEAhacks

## Overview
My mods to the [UCSF OPSCEA project](https://github.com/Kleen-Lab/OPSCEA) (March 2022 version).  Functionality added:
1. Using [Brainstorm](https://neuroimage.usc.edu/brainstorm/Introduction) output for electrode localization;
2. A function, `OPSCEA_struct`, that exports slice images (without EEG heatmaps, including both coronal and axial/sagittal views) and compiles a pdf document of the reconstruction; and
3. Simplified interface for OPSCEAparams.xls, that eliminates the need to specify the extreme contacts and desired color code for each electrode.

## Workflow
1. Run freesurfer.
2. In Brainstorm,
 * Load MRI and co-register CT;
 * Mark electrodes;
 * Export whole-brain figures of recon.
3. Run `create_opscea_subj.py`
4. Edit OPSCEAparams.xls
5. Run `OPSCEA.m` or `OPSCEA_struct.m`.

## Future Directions
1. Allow use of CAT12 recons;
2. Add [epileptogenicity score](https://pubmed.ncbi.nlm.nih.gov/18556663/) as alternate to line-length for EEG heatmaps;
3. Parallelize movie rendering.