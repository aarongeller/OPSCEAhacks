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

if ~exist('subj', 'var')
    % debug mode
    subj = 'notUCHZG';
    fssubjdir = fullfile(fsdir, 'UCHZG');
    eegdir = fullfile(fssubjdir, 'eeg');
else
    fssubjdir = fullfile(fsdir, subj);
    eegdir = fullfile(fssubjdir, 'eeg');
    if ~noeeg && ~norename
        % anonymize EDF file
        cd(eegdir);
        clean_edf_data(subj);
    end
end

ctdir = fullfile(fssubjdir, 'ct');
ctfile = dir(fullfile(ctdir, '*.nii.gz'));
eegfiles = dir(fullfile(eegdir, '*.edf'));

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
    for k=1:length(eegfiles)
        % review raw EEG file
        sFileRaw = bst_process('CallProcess', 'process_import_data_raw', [], [], ...
                               'subjectname',    subj, ...
                               'datafile',       {{fullfile(eegdir,eegfiles(k).name)}, 'SEEG-ALL'}, ...
                               'channelreplace', 0, 'channelalign',   0);

        % edit channelfile
        % 1) set type to SEEG
        bst_process('CallProcess', 'process_channel_setseeg', sFileRaw, []);

        % Get channel file
        [sStudy, iStudy] = bst_get('ChannelFile', sFileRaw.ChannelFile);
        ChannelMat = in_bst_channel(sFileRaw.ChannelFile);
        misc_channels = {'EKG1', 'Annotations', 'SpO2', 'EtCO2', 'Pulse', 'CO2Wave'};
        set_to_misc_str = '';
        isfirst = 1;
        for i=1:length(ChannelMat.Channel)
            % 2) find channels named EKG1, Annotations, SpO2, EtCO2, Pulse, CO2Wave
            for j=1:length(misc_channels)
                if strcmp(ChannelMat.Channel(i).Name, misc_channels(j))
                    if isfirst
                        isfirst = 0;
                    else
                        set_to_misc_str = [set_to_misc_str ', '];
                    end
                    set_to_misc_str = [set_to_misc_str ChannelMat.Channel(i).Name];
                    break;
                end
            end

            % 3) find channels starting with $, rename them and set to MISC
            if strcmp(ChannelMat.Channel(i).Name(1), "$")
                newname = ['XXX' int2str(i)];
                ChannelMat.Channel(i).Name = newname;
                if isfirst
                    isfirst = 0;
                else
                    set_to_misc_str = [set_to_misc_str ', '];
                end
                set_to_misc_str = [set_to_misc_str newname];
            end
        end

        bst_save(file_fullpath(sFileRaw.ChannelFile), ChannelMat, 'v7', 1);

        bst_process('CallProcess', 'process_channel_settype', sFileRaw, [], ...
                    'sensortypes', set_to_misc_str, 'newtype', 'MISC');
    end
end

toc;

% manual parts: 
% check coregistration
% mark electrodes
