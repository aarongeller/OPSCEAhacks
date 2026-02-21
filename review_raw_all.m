function review_raw_all(subj, norename)

fsdir = '/Applications/freesurfer/subjects';

fssubjdir = fullfile(fsdir, subj);
eegdir = fullfile(fssubjdir, 'eeg');

if ~exist('norename', 'var')
    norename = 0;
end

if ~norename
    % anonymize EDF file
    cd(eegdir);
    clean_edf_data(subj);
end

eegfiles = dir(fullfile(eegdir, '*.edf'));

for k=1:length(eegfiles)    
    % review raw EEG file
    review_raw_edf(subj, fullfile(eegdir, eegfiles(k).name));
end
