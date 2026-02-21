function create_brainstorm_subj(subj, noeeg, norename)

% prepare subject in brainstorm: anonymyze EDF, load MRI, coregister CT, load EEG
% and edit channel file

tic;

% check if brainstorm is running and if not, start it
if ~brainstorm('status')
    brainstorm;
end

if ~exist('noeeg', 'var')
    noeeg = 0;
end

if ~exist('norename', 'var')
    norename = 0;
end

fsdir = '/Applications/freesurfer/subjects';
fssubjdir = fullfile(fsdir, subj);
ctdir = fullfile(fssubjdir, 'ct');
ctfile = dir(fullfile(ctdir, '*.nii.gz'));

protocolname = 'IEEG_visualization';
gui_brainstorm('SetCurrentProtocol', bst_get('Protocol', protocolname));

% Create subject
[~, iSubject] = db_add_subject(subj, [], 0, 0);

% Import fs directory
isInteractive = 0;
nVertices = 15000;
errorMsg = import_anatomy_fs(iSubject, fssubjdir, nVertices, isInteractive, [], 0, 1);

% import CT and coregister
ProtocolSubjects = bst_get('ProtocolSubjects');
sSubject = ProtocolSubjects.Subject(iSubject);
refMriFile = sSubject.Anatomy(1).FileName;
postimplantct = import_mri(iSubject, fullfile(ctdir, ctfile.name), 'ALL', 0, 1);  
[postimplantctreg, errMsg, fileTag, sMriPostReg] = mri_coregister(postimplantct, refMriFile, 'spm', 1);
file_delete(file_fullpath(postimplantct), 1);
db_reload_subjects(iSubject);

if ~noeeg
    review_raw_all(subj, norename);
end

toc;

% manual parts: 
% check coregistration
% mark electrodes
