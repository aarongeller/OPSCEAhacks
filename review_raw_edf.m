function review_raw_edf(subj, fname)

% check if brainstorm is running and if not, start it
if ~brainstorm('status')
    brainstorm;
end

protocolname = 'IEEG_visualization';
gui_brainstorm('SetCurrentProtocol', bst_get('Protocol', protocolname));

% review raw EEG file
sFileRaw = bst_process('CallProcess', 'process_import_data_raw', [], [], ...
                       'subjectname',    subj, ...
                       'datafile',       {{fname}, 'SEEG-ALL'}, ...
                       'channelreplace', 0, 'channelalign',   0);

% edit channelfile

% Get channel file
[sStudy, iStudy] = bst_get('ChannelFile', sFileRaw.ChannelFile);
ChannelMat = in_bst_channel(sFileRaw.ChannelFile);
misc_channels = {'EKG', 'EKG1', 'Annotations', 'SpO2', ...
                 'EtCO2', 'Pulse', 'CO2Wave', 'BP1', 'BP2'};
set_to_misc_str = '';
misc_inds = [];
isfirst = 1;
for i=1:length(ChannelMat.Channel)
    % 1) find channels named EKG1, Annotations, SpO2, EtCO2, Pulse, CO2Wave
    for j=1:length(misc_channels)
        if strcmp(ChannelMat.Channel(i).Name, misc_channels{j})
            misc_inds(end+1) = i;
            if isfirst
                isfirst = 0;
            else
                set_to_misc_str = [set_to_misc_str ', '];
            end
            set_to_misc_str = [set_to_misc_str ChannelMat.Channel(i).Name];
            break;
        end
    end

    % 2) find channels starting with $, rename them and set to MISC
    if strcmp(ChannelMat.Channel(i).Name(1), "$")
        misc_inds(end+1) = i;
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

% 3) set channel Group
[ChannelMat.Channel(misc_inds).Group] = deal('MISC');
bst_save(file_fullpath(sFileRaw.ChannelFile), ChannelMat, 'v7');

% 4) set channel Type
bst_process('CallProcess', 'process_channel_settype', sFileRaw, [], ...
            'sensortypes', set_to_misc_str, 'newtype', 'MISC');
