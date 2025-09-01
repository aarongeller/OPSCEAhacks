function export_bs_figs(subj, channel_mat, outputdir, colors, labels)

% export figures from Brainstorm showing all electrodes, viewed
% from left, right, front and bottom, and export movies showing
% brain rotating around the Z axis and Y axis.

ch_mat_parts = split(channel_mat, '/');
protocolname = ch_mat_parts{end-4};

% check if brainstorm is running and if not, start it
if ~brainstorm('status')
    brainstorm nogui
end

% set electrode colors to match OPSCEA colors
load(channel_mat);
[~, labelsort] = sort(labels);
[~, iesort] = sort({IntraElectrodes.Name}); % usually will be sorted unless noedf
label_rows_to_skip = 3;
adjustnum = 0;
for i=1:length(IntraElectrodes)
    if i + label_rows_to_skip - adjustnum > length(labelsort)
        break;
    elseif ~strcmp(IntraElectrodes(iesort(i)).Name, ...
                   labels{labelsort(i + label_rows_to_skip - adjustnum)})
        % handle case where IntraElectrodes has non-SEEG channels
        adjustnum = adjustnum + 1;
    else
        IntraElectrodes(iesort(i)).Color = colors{labelsort(i + label_rows_to_skip - adjustnum)};
    end
end
s = who('-file', channel_mat);
save(channel_mat, s{:});

gui_brainstorm('SetCurrentProtocol', bst_get('Protocol', protocolname));
[hFig, iDS, iFig] = view_channels_3d({channel_mat}, 'SEEG', 'cortex', 1, 0);
set(hFig, 'Position', [10 10 692 572]);

views = {'left', 'right', 'bottom', 'front'};
output_prefix = 'SEEG_3D_';
movoptions.Duration = 10; % sec
movoptions.FrameRate = 15; % fps
movoptions.Quality = 100;

for i=1:length(views)
    figure_3d('SetStandardView', hFig, views{i});
    out_name = [outputdir '/' output_prefix subj '_' views{i}];
    out_figure_image(hFig, [out_name '.png']);
    if strcmp(views{i}, 'bottom') || strcmp(views{i}, 'front')
        out_figure_movie(hFig, [out_name '.avi'], 'horizontal', movoptions);
    end
end

close all;
