function setup_bs(subj)

% prepare subject in brainstorm: load MRI, coregister CT, load EEG
% and edit channel file

% check if brainstorm is running and if not, start it
if ~brainstorm('status')
    brainstorm;
end

fsdir = '/Applications/freesurfer/subjects';

if ~exist('subj', 'var')
    subj = 'notUCHZG';
    fssubjdir = fullfile(fsdir, 'UCHZG');
else
    fssubjdir = fullfile(fsdir, subj);
end

ctdir = fullfile(fssubjdir, 'ct');
ctfile = dir(fullfile(ctdir, '*.nii.gz'));
eegdir = fullfile(fssubjdir, 'eeg');
eegfiles = dir(fullfile(eegdir, '*.edf'));

protocolname = 'IEEG_Visualization';
gui_brainstorm('SetCurrentProtocol', bst_get('Protocol', protocolname));

% % anonymize EDF file
% cd eegdir;
% clean_edf_data(subj);

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

for k=1:length(eegfiles)
    % review raw EEG file
    sFileRaw = bst_process('CallProcess', 'process_import_data_raw', [], [], ...
                           'subjectname',    subj, ...
                           'datafile',       {{fullfile(eegdir,eegfiles(k).name)}, 'SEEG-ALL'}, ...
                           'channelreplace', 0, ...
                           'channelalign',   0);

    % edit channelfile
    % 1) set type to SEEG
    bst_process('CallProcess', 'process_channel_setseeg', sFileRaw, []);
    % 2) find channels named EKG1, Annotations, SpO2, EtCO2, Pulse, CO2Wave
    % Get channel file
    [sStudy, iStudy] = bst_get('ChannelFile', sFileRaw.ChannelFile);
    ChannelMat = in_bst_channel(sFileRaw.ChannelFile);
    misc_channels = {'EKG1', 'Annotations', 'SpO2', 'EtCO2', 'Pulse', ...
                     'CO2Wave'};
    misc_present_str = '';
    isfirst = 1;
    for i=1:length(ChannelMat.Channel)
        for j=1:length(misc_channels)
            if strcmp(ChannelMat.Channel(i).Name, misc_channels(j))
                if isfirst
                    isfirst = 0;
                else
                    misc_present_str = [misc_present_str ', '];
                end
                misc_present_str = [misc_present_str ChannelMat.Channel(i).Name];
                break;
            end
        end
    end

    % 3) find channels starting with $, rename them and set to MISC

    dollar_chan_str = '';
    isfirst = 1;
    for i=1:length(ChannelMat.Channel)
        for j=1:length(misc_channels)
            if strcmp(ChannelMat.Channel(i).Name(1), "$")
                newname = ['XXX' int2str(i)];
                ChannelMat.Channel(i).Name = newname;
                if isfirst
                    isfirst = 0;
                else
                    dollar_chan_str = [dollar_chan_str ', '];
                end
                dollar_chan_str = [dollar_chan_str newname];
                break;
            end
        end
    end

    bst_save(file_fullpath(sFileRaw.ChannelFile), ChannelMat, 'v7', 1);

    bst_process('CallProcess', 'process_channel_settype', sFileRaw, [], ...
                'sensortypes', misc_present_str, ...
                'newtype',     'MISC');

    bst_process('CallProcess', 'process_channel_settype', sFileRaw, [], ...
                'sensortypes', dollar_chan_str, ...
                'newtype',     'MISC');
end

% manual parts: 
% check coregistration
% mark electrodes
