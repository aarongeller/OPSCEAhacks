function OPSCEAsurfslice_axial(subject,orientation,elecs,weights,subj_dir,fs_dir,cax,CM,gsp,j,isfirstframe,adjust_coords,force_angle)
%     (This is a subfunction created as a part of) Omni-planar and surface
%     casting of epileptiform activity (OPSCEA) (UC Case Number SF2020-281)
%     jointly created by Dr. Jon Kleen, Ben Speidel, Dr. Robert Knowlton,
%     and Dr. Edward Chang is licensed for non-commercial research use at
%     no cost by the Regents of the University of California under CC
%     BY-NC-SA 4.0 (https://creativecommons.org/licenses/by-nc-sa/4.0/).
%     Please contact innovation@ucsf.edu if you are interested in using
%     OPSCEA for commercial purposes.
%     The following copyright notice and citation is to be included in any
%     publication, material or media wherein all or a part of Licensed
%     Material is contained, “Certain materials incorporated herein are
%     Copyright © 2016 The Regents of the University of California
%     (REGENTS). All Rights Reserved.
%     Please cite the following paper in your publications if you have used
%     our software in your research, as well as any relevant toolboxes used
%     herein as appropriate (img_pipe, FreeSurfer): Kleen JK, Speidel B,
%     Baud MO, Rao VR, Ammanuel SG, Hamilton LS, Chang EF, Knowlton RC.
%     Accuracy of omni-planar and surface casting of epileptiform activity
%     for intracranial seizure localization. In press at Epilepsia.”
% Global variables, added to improve computation time (most need to be calculated the first frame only
global loaf; %contains data & info about the volume and surfaces
global sliceinfo; %contains data & info about each individual slice (for each depth)
global S; %contains data & info about the ECoG data and add'l parameters from excel sheet
global I; %contains elecmatrix (coordinates), weights (for heatmap), and index of channels to skip
global salphamask;

if nargin<3; error('surfslice requires at least 3 input arguments'); end

if ~exist('adjust_coords', 'var')
    adjust_coords = [0 0 0];
end

if ~exist('force_angle', 'var')
    force_angle = NaN;
end

%% this section sets up paths and filenames
fsbin = strcat(fs_dir,'/bin/');
subjectpath = strcat(subj_dir,'/',subject,'/Imaging');
brainpathmgz = strcat([subjectpath,'/mri/brain.mgz']);
aparcpathmgz = strcat([subjectpath,'/mri/aparc+aseg.mgz']);
bashpath=getenv('PATH');
newpath=strcat([bashpath, ':',fsbin]); 
if length(newpath)>1000; newpath='/Applications/freesurfer/bin/'; end
%%  this section loads in and sets up volume and surface data
setenv('BRAINMGZ',brainpathmgz);
if isempty(loaf.vrf)
    v=load_mgh(brainpathmgz); 
    apas = load_mgh(aparcpathmgz); 
    vr=imrotate3(v,90,[1 0 0],'cubic','crop');
    vr=imrotate3(vr,90,[0 1 0],'cubic','crop');
    vrf=flip(vr,2);
    vrf=flip(vrf,3);
    loaf.vrf=vrf;
else 
    vrf=loaf.vrf;
end
if isempty(loaf.apasrf) 
    apasr=imrotate3(apas,90,[1 0 0],'nearest','crop');
    apasr=imrotate3(apasr,90,[0 1 0],'nearest','crop');
    apasrf = flip(apasr,2);
    apasrf = flip(apasrf,3);
    loaf.apasrf=apasrf; clear apasrf
end
alphamask = ones(size(loaf.apasrf)); % *note: consider making alpha mask also clip points outside of axislim boundaries
alphamask(loaf.apasrf == 15|loaf.apasrf == 46|loaf.apasrf ==7|loaf.apasrf ==16|loaf.apasrf == 8|loaf.apasrf == 47|loaf.apasrf == 0)=0;
hold on;

if isfirstframe
    [m,b,theta,e1] = get_mb(elecs, [1 0 1], force_angle); % get line in XZ plane
    intrcpt = b + 128;
    maxgrad = get_max_gradient(elecs);
    if maxgrad==3
        orientation = 's';
        sliceinfo(j).sagittal = 1;
    else
        sliceinfo(j).sagittal = 0;
    end
    YY = meshgrid(-128:128)'; % coordinates for spatial calculations and plotting
    yslice = meshgrid(1:256)' + adjust_coords(2);
    xslice = cos(theta).*meshgrid(-127.5:127.5) + elecs(e1,1) + 128.5 + adjust_coords(1);
    zslice = sin(theta).*meshgrid(-127.5:127.5) + elecs(e1,3) + 128.5 + adjust_coords(3);
    XX = meshgrid(-128:128).*cos(theta) + elecs(e1,1);
    ZZ = sin(theta).*meshgrid(-128:128) + elecs(e1,3); 
    sliceinfo(j).azel = [180, -90]; % for pure inferior view: [180, -90]
    sliceinfo(j).XX = XX; 
    sliceinfo(j).YY = YY; 
    sliceinfo(j).ZZ = ZZ; 
    sliceinfo(j).xslice = xslice; 
    sliceinfo(j).yslice = yslice; 
    sliceinfo(j).zslice = zslice;
    sliceinfo(j).viewangle = rad2deg(theta);
    sliceinfo(j).b = b;
    % sliceinfo(j).sl=m; 
    sliceinfo(j).slicenum=intrcpt; 
    %create surface meshes that are split along the slice plane and plot
    %only one segment so that the sliceplane is still visible
    sliceinfo(j).lsplit = splitbrain(loaf.lpial, orientation, b, m);
    sliceinfo(j).rsplit = splitbrain(loaf.rpial, orientation, b, m);
    sliceinfo(j).final_orientation = orientation;
end
hold on;
if ~any(weights(:))
    % doing structural plot without heatmap
    if ~isempty(sliceinfo(j).lsplit.vert)
        % can happen for sagittal cuts that exclude one hemisphere
        lbrn=ctmr_gauss_plot_edited(sliceinfo(j).lsplit,[],[],S.cax,0,S.cm,S.gsp);
    end
    if ~isempty(sliceinfo(j).rsplit.vert)
        % can happen for sagittal cuts that exclude one hemisphere
        rbrn=ctmr_gauss_plot_edited(sliceinfo(j).rsplit,[],[],S.cax,0,S.cm,S.gsp);
    end
else
    lbrn=ctmr_gauss_plot_edited(sliceinfo(j).lsplit,I.em(I.nns,:),I.w8s(I.nns),S.cax,0,S.cm,S.gsp); 
    rbrn=ctmr_gauss_plot_edited(sliceinfo(j).rsplit,I.em(I.nns,:),I.w8s(I.nns),S.cax,0,S.cm,S.gsp); 
end
if isfirstframe
    s=slice(vrf,sliceinfo(j).xslice,sliceinfo(j).yslice,sliceinfo(j).zslice); 
    s.Visible = 'off';
    sliceinfo(j).CData = s.CData;
end
bread=double(sliceinfo(j).CData); % the slice
%This method tries to keep all slices the same brightness and contrast levels
if isempty(loaf.normloaf) 
    loaf.normloaf=prctile(reshape(vrf,1,prod(size(vrf))),99.9); %99.9th percentile, since outliers can give poor max
end
%improve contrast/brightness for plotting, %otherwise can be too dark and will obscure heatmap for some patients
bread=round(bread/loaf.normloaf*254); 
cim = toaster(sliceinfo(j).XX,sliceinfo(j).YY,sliceinfo(j).ZZ,bread,elecs,weights,[0 cax(2)],CM,gsp);
INPUT = uint8(cim); 
sliceimage = padarray(INPUT,[1 1],0,'post'); %padding to appropriate size
%Use the same slicing parameters to create an alphamask to make slice
%transparent outside of brainvolume
a = slice(alphamask,sliceinfo(j).xslice, sliceinfo(j).yslice, sliceinfo(j).zslice, 'nearest');
a.Visible = 'off';
AA = padarray(a.CData, [1 1],0, 'post');
AAnonan=AA; AAnonan(isnan(AA))=0; 
SE = strel('disk',2);
alphamap = bwareaopen(imopen(AAnonan,SE),50); 
[xedge, yedge, zedge] = getEdges(alphamap, sliceinfo(j).XX, sliceinfo(j).YY, sliceinfo(j).ZZ, 'a');
sliceinfo(j).corners=[xedge fliplr(xedge); yedge([1 1 2 2]); zedge fliplr(zedge);];   %for oblique slice planes
%create surface in coordinate space that slices brain in the
%appropriate plane and apply color and transparency data
surface(sliceinfo(j).XX,sliceinfo(j).YY,sliceinfo(j).ZZ,'CData',sliceimage,'EdgeColor','none','FaceColor','texturemap','FaceAlpha','texturemap','EdgeAlpha',0,'AlphaData',alphamap,'specularexponent',5);
shading flat
caxis(S.cax) % colormap(S.cm); 
alim([0.1 1])
set(gca,'Clipping','off')
axis vis3d
salphamask = alphamask;
