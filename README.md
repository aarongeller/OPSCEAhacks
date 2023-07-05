# OPSCEAhacks

## Overview
My mods to the [UCSF OPSCEA project](https://github.com/Kleen-Lab/OPSCEA) (March 2022 version).  

### Main functionality added:
1. Using [Brainstorm](https://neuroimage.usc.edu/brainstorm/Introduction) output for electrode localization;
2. A function, `OPSCEA_recon`, that exports slice images (without EEG heatmaps, including both coronal and axial/sagittal views) and compiles a pdf document of the reconstruction.

### Minor tweaks:
1. Simplified interface for OPSCEAparams.xls, that eliminates the need to specify the extreme contacts and desired color code for each electrode.
2. Movie filename includes seizure time (taken from seizure .mat name) and fps, and if an existing movie with the designated name already exists, a new name will be used to prevent automatic over-writing.

## Dependencies
* [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/)
* [MATLAB](https://www.mathworks.com/products/matlab.html) (with [Brainstorm](https://neuroimage.usc.edu/brainstorm/Introduction) and [OPSCEA](https://github.com/Kleen-Lab/OPSCEA))
* [Python](https://www.python.org/) (with [NiBabel](https://nipy.org/nibabel/) and [MNE](https://mne.tools/stable/index.html))
* [LaTeX](https://www.python.org/)

## Workflow
1. Run freesurfer.
2. In Brainstorm,
 * Load MRI,
 * Co-register CT, and
 * Mark electrodes.
3. Run `create_opscea_subj.py`
4. Edit OPSCEAparams.xls
5. Run `OPSCEA.m` or `OPSCEA_recon.m`.

## Future Directions
1. Allow use of [CAT12](https://neuro-jena.github.io/cat/) recons;
2. Add [epileptogenicity score](https://pubmed.ncbi.nlm.nih.gov/18556663/) as alternate to line-length for EEG heatmaps;
3. Parallelize movie rendering.
