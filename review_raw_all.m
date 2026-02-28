function review_raw_all(subj, protocolname, norename)

if ~exist('protocolname', 'var')
    protocolname = 'IEEG_visualization';
end

if ~exist('norename', 'var')
    norename = 0;
end

fsdir = '/Applications/freesurfer/subjects';

fssubjdir = fullfile(fsdir, subj);
eegdir = fullfile(fssubjdir, 'eeg');
bsdatadir = fullfile('/Users/aaron/Documents/brainstorm_db/', ...
                     protocolname, 'data');
subjbsdatadir = fullfile(bsdatadir, subj);

if ~norename
    % anonymize EDF file
    cd(eegdir);
    clean_edf_data(subj);
end

eegfiles = dir(fullfile(eegdir, '*.edf'));

for k=1:length(eegfiles)
    thisfile = eegfiles(k).name;
    [~, fstem, ~] = fileparts(thisfile);

    % check that it's not already imported
    if exist(fullfile(subjbsdatadir, ['@raw' fstem]), 'dir')
        display(['Skipping ' thisfile '...']);
    else
        % review raw EEG file
        review_raw_edf(subj, fullfile(eegdir, thisfile));
    end
end
