function export_bs_figs(subj)

% export figures from Brainstorm showing all electrodes, viewed
% from left, right, front and bottom, and export movies showing
% brain rotating around the Z axis and Y axis.

fs_output_path = ['/Applications/freesurfer/subjects/' subj '/ct'];
protocolname = 'IEEG_visualization';
bs_subj_path = ['/Users/aaron/Documents/brainstorm_db/' protocolname '/data/' subj];
chfiles = dir([bs_subj_path '/@raw*']);
FileName = [subj '/' chfiles(1).name '/channel.mat'];

gui_brainstorm('SetCurrentProtocol', bst_get('Protocol', protocolname));
[hFig, iDS, iFig] = view_channels_3d({FileName}, 'SEEG', 'cortex', 1, 0);

views = {'left', 'bottom', 'right', 'front'};
output_prefix = 'SEEG_3D_';
movoptions.Duration = 10; % sec
movoptions.FrameRate = 15; % fps
movoptions.Quality = 75;

for i=1:length(views)
    figure_3d('SetStandardView', hFig, views{i});
    out_name = [fs_output_path '/' output_prefix subj '_' views{i}];
    out_figure_image(hFig, [out_name '.png']);
    if strcmp(views{i}, 'bottom') || strcmp(views{i}, 'front')
        out_figure_movie(hFig, [out_name '.avi'], 'horizontal', movoptions);
    end
end

close all;
