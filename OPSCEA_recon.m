function OPSCEA_recon(pt, selected_leads, force_angle_coronal, force_angle_axial, dopdf)
% Do a limited (structure only) OPSCEA run for SEEG electrodes,
% exporting coronal and axial views together with a corresponding
% surface view illustrating the cutplane, then calling
% make_slice_pdf.py to collect them in a pdf.
%
% EXAMPLE USAGE: OPSCEA_struct(pt, dopdf, showlabels)
% 
% pt is a string such as 'UCSF4' or 'JaneDoe', acts as a prefix for files below
% dopdf is 1 (default) or 0
% showlabels is:
%      1 if you want to show the channel labels (default)
%      0 to hide them AND randomize channels (blinding the reader to the
%      electrode locations of the trace-based ICEEG as in Kleen et al. 2021)

%     Omni-planar and surface casting of epileptiform activity (OPSCEA) (UC
%     Case Number SF2020-281) jointly created by Dr. Jon Kleen, Ben
%     Speidel, Dr. Robert Knowlton, and Dr. Edward Chang is licensed for
%     non-commercial research use at no cost by the Regents of the
%     University of California under CC BY-NC-SA 4.0
%     (https://creativecommons.org/licenses/by-nc-sa/4.0/). Please contact
%     innovation@ucsf.edu if you are interested in using OPSCEA for
%     commercial purposes.

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

% needs: statistics toolbox, signal processing toolbox, image
% processing toolbox

if ~exist('showlabels','var')||isempty(showlabels)
    % default displays ICEEG and depth labels
    showlabels=true; 
end 

if ~exist('dopdf', 'var')
    dopdf = 1;
end

sz = '01'; % remove

opsceapath=['/Users/aaron/Documents/MATLAB/OPSCEA-main/']; % AG

% opsceapath=['/Users/kleentestaccount/Desktop/OPSCEA/'];   %path for parameters sheet
opsceadatapath=[opsceapath 'OPSCEADATA/'];   %path for OPSCEA ICEEG and imaging data
if ~exist(opsceadatapath,'dir')
    error('Directory for your data needs to be corrected'); 
end
cd(opsceapath);

ptsz=[pt '_' sz]; % prefix for filenames of specific seizure
ptpath=[opsceadatapath pt '/']; % patient's folder
szpath= [ptpath ptsz '/']; % specific seizure's folder
disp(['Running ' pt ', seizure ' sz '...']);

%% Initiate global variables
  global S; % holds general parameters
 % for speed, these are filled during first frame of surfslice then re-used
  global loaf; 
  global sliceinfo; 
  global I;


%% Import parameters
% for specific seizure 
[~,prm_allPtSz]=xlsread([opsceapath 'OPSCEAparams'],'params'); 
fields_SZ=prm_allPtSz(1,:); % header for columns of seizure parameters
prm=prm_allPtSz(strcmp(pt,prm_allPtSz(:,1))&strcmp(sz,prm_allPtSz(:,2)),:);
if isempty(prm); error(['ATTENTION: No entry exists for ' pt ' seizure ' sz ' in the params master sheet']); end
% Import parameters for patient's specific plot (layout of video frame)
[~,plt]=xlsread([opsceapath 'OPSCEAparams'],pt); 
fields_PLOT=plt(1,:); plt(1,:)=[]; % header for columns of plotting parameters
plottype=plt(:,strcmpi(fields_PLOT,'plottype')); %type of plot for each subplot (accepts: iceeg, surface, depth, or colorbar)

if ~exist('selected_leads', 'var')
    rows_to_do = 1:size(plt,1);
    force_angle_coronal = nan(1, size(plt,1));
    force_angle_axial = nan(1, size(plt,1));
else
    rows_to_do = get_matching_labels(plt(:,end), selected_leads);
    num_to_do = length(rows_to_do) - 3;

    if exist('force_angle_coronal', 'var')
        if length(force_angle_coronal)>0
            if length(force_angle_coronal)~=num_to_do
                error(['Angle list lengths must match electrode list lengths.  Quitting.']);
            end
            nanvec = nan(1, max(rows_to_do));
            for i=1:num_to_do
                nanvec(rows_to_do(i+3)) = force_angle_coronal(i);
            end
            force_angle_coronal = nanvec;
        else
            force_angle_coronal = nan(1, max(rows_to_do));
        end
    else
        force_angle_coronal = nan(1, max(rows_to_do));
    end

    if exist('force_angle_axial', 'var')
        if length(force_angle_axial)>0
            if length(force_angle_axial)~=num_to_do
                error(['Angle list lengths must match electrode list lengths.  Quitting.']);
            end
            nanvec = nan(1, max(rows_to_do));
            for i=1:num_to_do
                nanvec(rows_to_do(i+3)) = force_angle_axial(i);
            end
            force_angle_axial = nanvec;
        else
            force_angle_axial = nan(1, max(rows_to_do));
        end
    else
        force_angle_axial = nan(1, max(rows_to_do));
    end
end

cd 
%% prepare subplot specifications
subplotrow=str2double(plt(:,strcmpi(fields_PLOT,'subplotrow')));
subplotcolumn=str2double(plt(:,strcmpi(fields_PLOT,'subplotcolumn')));
subplotstart=plt(:,strcmpi(fields_PLOT,'subplotstart')); 
subplotstop=plt(:,strcmpi(fields_PLOT,'subplotstop')); 
for j=1:length(plottype); 
    subplotnum{j,1}=str2double(subplotstart{j}):str2double(subplotstop{j});
end
surfaces=plt(:,strcmpi(fields_PLOT,'surfaces'));
surfacesopacity=plt(:,strcmpi(fields_PLOT,'surfacesopacity'));
viewangle=lower(plt(:,strcmpi(fields_PLOT,'viewangle')));

%% parcel all individual depth labels, contact #s, and colors. If no depths, make it  =[];
depthlabels=plt(:,strcmpi(fields_PLOT,'depthlabels'));
isdepth=strcmpi(plottype,'depth'); 
depths=cell(size(isdepth));
if any(isdepth)
    % check for contactdict.mat
    if exist([ptpath '/Imaging/Elecs/contactdict.mat'], 'file')
        % if present, load and convert to dict
        load([ptpath '/Imaging/Elecs/contactdict.mat']);
        contactmap = containers.Map();
        for i=1:length(wirelabels)
            contactmap(wirelabels{i}) = extremecontacts(i,:);
        end
        depths = {};
        % make depths array
        numnans = 0;
        for i=1:length(depthlabels)
            if isempty(depthlabels{i})
                depths{i} = nan;
                depthcolor{i} = nan;
                numnans = numnans + 1;
            else
                thiselec = contactmap(depthlabels{i});
                depths{i} = thiselec(1):thiselec(2);           
                depthcolor{i} = colors(i-numnans,:);
            end
        end
        depths = depths';
        depthcolor = depthcolor';

    else
        depthEfirst=plt(:,strcmpi(fields_PLOT,'depthEfirst')); 
        depthElast=plt(:,strcmpi(fields_PLOT,'depthElast')); 
        for j=1:length(depths); 
            depths{j}=str2double(depthEfirst{j}):str2double(depthElast{j});
        end
        depthcolor=plt(:,strcmpi(fields_PLOT,'depthcolor')); 
        for j=1:length(depthcolor)
            splt=regexp(depthcolor{j},',','split'); 
            depthcolor{j}=str2double(splt); 
        end         
    end
    pltzoom=str2double(plt(:,strcmpi(fields_PLOT,'pltzoom')));
    pltshowplanes=str2double(plt(:,strcmpi(fields_PLOT,'showplanes')))==1; %logical index of plots in which to show slice planes
end

%% Get time segments within the ICEEG file to use
VIDstart=prm(:,strcmpi(fields_SZ,'VIDstart')); VIDstop=prm(:,strcmpi(fields_SZ,'VIDstop')); %chunk of data (seconds into ICEEG data file) to use from the whole ICEEG data clip for the video
S.VIDperiod=[str2double(VIDstart{1}) str2double(VIDstop{1})];
if exist('timewindow', 'var') & ~isempty(timewindow)
    S.VIDperiod = timewindow;
end
BLstart=prm(:,strcmpi(fields_SZ,'BLstart')); BLstop=prm(:,strcmpi(fields_SZ,'BLstop')); %chunk of data (seconds into ICEEG data file) to use for baseline (for z-score step)
S.BLperiod=[str2double(BLstart{1}) str2double(BLstop{1})];

%transform, scaling, and display options
S.llw=str2double(prm{strcmp('llw',fields_SZ)}); %default linelength window (in seconds)
S.iceeg_scale=prm{strcmp('iceeg_scale',fields_SZ)}; %percentile (number >50 and <100), used here similar to gain ICEEG waveform display, usually 95
if ischar(S.iceeg_scale); S.iceeg_scale=str2double(S.iceeg_scale); end 
S.fps=str2double(prm{strcmp('fps',fields_SZ)});             %frames per sec of ICEEG (default 15)
S.cax=str2double(regexp(prm{strcmp('cax',fields_SZ)},',','split'));         %color axis for heatmap
S.gsp=str2double(prm{strcmp('gsp',fields_SZ)}); %gaussian spreading parameter (default 10)
params={'iceeg_scale','fps','cax','gsp'}; 
paramsnans=isnan([(isnan(S.iceeg_scale) | S.iceeg_scale<=50 | S.iceeg_scale>=100)   S.fps   any(isnan(S.cax)) S.gsp]); 
if any(paramsnans); error(['ATTENTION OPSCEA USER: The "' params{paramsnans} '" term(s) is/are in an incorrect format (perhaps number instead of string), check excel seizure parameter sheet']); 
end 
cm=prm{strcmp('cm',fields_SZ)};
switch cm; case 'cmOPSCEAcool'; cm=cmOPSCEAcool; 
  case 'cmOPSCEAjet'; cm=cmOPSCEAjet; 
end
S.cm=cm; %colormap to use for heatmap
S.iceegwin=str2double(prm{strcmp('iceegwin',fields_SZ)}); %how much trace-based ICEEG to view at a time in the ICEEG window
S.marg=str2double(prm{strcmp('marg',fields_SZ)}); %offset of real-time LL txform from beginning of viewing window (in sec; converts to samples below)
S.slicebright=str2double(prm{strcmp('slicebright',fields_SZ)}); if isnan(S.slicebright); S.slicebright=0; end %brighten up slices (usually 0 to 50)


% additional adjustment for display window
S.VIDperiod=[S.VIDperiod(1)-S.marg   S.VIDperiod(2)+S.iceegwin-S.marg]; 
S.fields=fields_SZ; clear fields

S.prm=prm; clear prm
S.prm_allPtSz=prm_allPtSz; clear prm_allPtSz
S=orderfields(S); %alphabetize the structure fields for ease of use/search

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% load ICEEG data, and the bad channels verified for that specific data
szfname = dir([szpath ptsz '*_sz.mat']).name;
load([szpath szfname])
% load([szpath ptsz '_badch']); 
% if size(d,1)>size(d,2); d=d'; end % orient to channels by samples
[nch,ntp]=size(d); f=1; 
% disp(['Length of data to play for video is ' num2str(round(ntp/sfx)) ' sec'])

% % error checks for selected time periods
% if any([S.VIDperiod(1) S.BLperiod(1)]<0)
%     error('VIDperiod is out of bounds of file (time < 0). Check both VIDstart and BLstart times and make sure the "marg" value (subtracted from VIDstart and BLstart), cannot be < 0'); 
% elseif any([S.VIDperiod(2) S.BLperiod(2)]>ntp)
%     error('VIDperiod is beyond the length of the file. Check VIDstop and BLstop times vs. length of actual ICEEG data'); 
% end 

%% locate and load electrode file for labels and XYZ coordinates
load([ptpath 'Imaging/Elecs/Electrodefile.mat']); 
if ~exist('anatomy','var'); anatomy=cell(size(elecmatrix,1),4); end
if size(anatomy,1)>size(elecmatrix,1); anatomy(size(elecmatrix,1)+1:end)=[]; end
anat=anatomy; clear anatomy; if size(anat,2)>size(anat,1); anat=anat'; end
if size(anat,2)==1; anat(:,2)=anat(:,1); end; 
if ~exist('eleclabels','var'); eleclabels=anat(:,1); end
em=elecmatrix; clear elecmatrix; emnan=isnan(mean(em,2)); badch(emnan)=1; em(emnan,:)=0; EKGorREF=strcmpi('EKG1',anat(:,1))|strcmpi('EKG2',anat(:,1))|strcmpi('EKG',anat(:,2))|strcmpi('EKGL',anat(:,2))|strcmpi('REF',anat(:,1)); anat(EKGorREF,:)=[]; em(EKGorREF,:)=[]; eleclabels(EKGorREF,:)=[]; 
   
%% load meshes you want to plot
meshpath='Imaging/Meshes/';
Rcortex=load([ptpath meshpath pt '_rh_pial.mat']); 
loaf.rpial=Rcortex; 
Rcrtx=Rcortex.cortex; 
clear Rcortex

Lcortex=load([ptpath meshpath pt '_lh_pial.mat']); 
loaf.lpial=Lcortex; 
Lcrtx=Lcortex.cortex; 
clear Lcortex

for i=1:length(surfaces) 
    hippentry(i)=~isempty(strfind(surfaces{i},'hipp')); 
    amygentry(i)=~isempty(strfind(surfaces{i},'amyg')); 
end 
errmsg='ATTN: MISSING A MESH, need to add this mesh file to directory (or remove/omit from frame): ';

if any(hippentry)
    Rhipp=[ptpath meshpath 'subcortical/rHipp_subcort.mat']; 
    Lhipp=Rhipp; 
    Lhipp(end-16)='l'; 
    if exist(Rhipp,'file')
        Rhipp=load(Rhipp); 
        Rhipp=Rhipp.cortex;  
        Lhipp=load(Lhipp); 
        Lhipp=Lhipp.cortex;     
    else 
        error([errmsg 'hipp']); 
    end 
end

if any(amygentry)
    Ramyg=[ptpath meshpath 'subcortical/rAmgd_subcort.mat']; 
    Lamyg=Ramyg; 
    Lamyg(end-16)='l'; 
    if exist(Ramyg,'file')
        Ramyg=load(Ramyg); 
        Ramyg=Ramyg.cortex;  
        Lamyg=load(Lamyg); 
        Lamyg=Lamyg.cortex;     
    else 
        error([errmsg 'amyg']); 
    end 
end

if exist('Lhipp', 'var') && exist('Lamyg', 'var') && exist('extendedObjectMesh')
    Lcrtx = undo_fs_ahectomy(Lcrtx, Lhipp, Lamyg);
    loaf.lpial.cortex = Lcrtx;
end

if exist('Rhipp', 'var') && exist('Ramyg', 'var') && exist('extendedObjectMesh')
    Rcrtx = undo_fs_ahectomy(Rcrtx, Rhipp, Ramyg);
    loaf.rpial.cortex = Rcrtx;
end

drows=find(strcmp(plottype,'depth'))'; ndepths=length(drows);
depthch=[]; 
for i=1:length(drows) 
    depthch=[depthch depths{drows(i)}]; 
end; 
clear i %identify all depth electrode channels

isR=em(:,1)>0; 
isR=nansum(em(:,1))>0; % doesn't make sense to take sum over all contacts...
isL=isR~=1; %handy binary indicators for laterality

isRdepth = [];
isLdepth = [];
for i=1:length(depths)
    if ~isnan(depths{i})        
        xval_highcontact = em(depths{i}(end),1);
        isRdepth(end+1) = xval_highcontact>=0;
        isLdepth(end+1) = xval_highcontact<0;
    else
        isRdepth(end+1) = nan;
        isLdepth(end+1) = nan;
    end
end

%% get xyz limits for plotting purposes
perim=1; % how many millimeters away from brain/electrodes boundaries to set the colorcoded plane perimeter, recommend >0 to avoid skimming brain surface (default 1mm)
axl(:,:,1)=[min([Rcrtx.vert; Lcrtx.vert]); max([Rcrtx.vert; Lcrtx.vert])]; %min and max of MESH VERTICES' x y z coordinates (2x3)
axl(:,:,2)=[min(em); max(em)]; %%min and max of ELECTRODES' x y z coordinates (2x3) and stack them (2x3x2)
axl=[min(axl(1,:,:),[],3); max(axl(2,:,:),[],3)]; % Get the minima and maxima of both
axislim=reshape(axl,1,6)+[-1 1 -1 1 -1 1]*perim; clear axl %Use the, to define the axis boundaries, and add additional perimeter (perim)

% %% formatting checks, and consolidation of bad channels
ns=unique( [find(badch);   find(isnan(mean(em,2)));   find(isnan(mean(d,2)))]  ); % bad channels: those that are pre-marked, or if NaNs in their coordinates or ICEEG data traces
nns=true(nch,1); nns(ns)=0; %nns=find(nns); %consolidate bad channels and those with NaNs
% %remove data channels previously clipped at end. Only include that which has electrode coordinates (intracranial)
% if size(em,1)>size(d,1); nch=size(em,1); d(nch+1:end,:)=[]; LL(nch+1:end,:)=[]; nns(nch+1:end)=[]; ns(ns>size(em,1))=[]; 
%    fprintf(2, 'ALERT: Clipping off extra bad channel entries (make sure you have the right ICEEG and bad channel files loaded)\n');
% end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTING TIME!
sliceinfo=[]; loaf.vrf=[]; loaf.apasrf=[]; loaf.normloaf=[]; ...
sliceinfo.viewangle=zeros(size(plt,1),3); sliceinfo.azel=[]; ...
sliceinfo.corners=[]; loaf.isR=isR; loaf.isL=isL; loaf.isRdepth=isRdepth; loaf.isLdepth=isLdepth;
clear F; 
nch=length(find(nns)); 

chanorder=1:size(d(nns,:),1); 

figure(1);
set(gcf, 'color','w');
set(gcf, 'Position',[1 5 1280 700]);

planes = {'a', 'c'};
w8s = [];

for i=1:length(planes)
    S.sliceplane = planes{i};
    subplot(1,1,1); %clears all axes, to start fresh each frame
    for j=rows_to_do
        h = subplot(subplotrow(j),subplotcolumn(j),subplotnum{j}); 
        switch upper(plottype{j,1})
          case 'SURFACE' % plotting surfaces only
            hold off; srf=regexp(surfaces{j},',','split'); % list the specific surfaces wanted for this subplot
            srfalpha=regexp(surfacesopacity{j},',','split'); % list their corresponding opacities (values from 0 to 1; 0=invisible, 1=opaque)
            if length(srf)~=length(srfalpha)
                ME = MException("BadParams:mismatch", 'Number of surface to plot does not match number of alpha designations, check excel sheet');
                throw(ME);
            end
            acceptedterms={'rcortex','lcortex','rhipp','lhipp','ramyg','lamyg','wholebrain'};
            for s=1:length(srf)
                srf{s}=lower(srf{s}); %convert to lower case for easier string matching
                if ~isempty(intersect(srf{s},acceptedterms)) %make sure user specified accepted terminologies for the meshes
                    switch srf{s}; %see below for case "wholebrain"s
                      case 'rcortex'; srfplot=Rcrtx; 
                      case 'lcortex'; srfplot=Lcrtx; 
                      case 'rhipp';   srfplot=Rhipp; 
                      case 'lhipp';   srfplot=Lhipp; 
                      case 'ramyg';   srfplot=Ramyg; 
                      case 'lamyg';   srfplot=Lamyg; 
                    end
                else
                    disp(['ATTN: Row ' num2str(j) ' defined as surface but does not contain an accepted mesh term']); 
                    disp(acceptedterms); 
                    error('');
                end
                % plot the individual heatmapped surface
                if exist('srfplot','var')
                    hh=ctmr_gauss_plot_edited(srfplot, [], [], S.cax, 0, S.cm, S.gsp); 
                    alpha(hh,str2double(srfalpha{s})); % Adjust opacity specified for that row
                else
                    disp(['ALERT: One of the entries in row ' num2str(j) ' is not a valid entry, accepts:']); 
                    disp(acceptedterms); 
                end
            end
            if isempty(intersect(srf{s},{'rcortex','lcortex'}))||strcmpi(srf,'wholebrain') %for glass brain (hipp and/or amyg only) and wholebrain plots
                glass1 = ctmr_gauss_plot_edited(Rcrtx, [], [], S.cax, 0, S.cm, S.gsp); 
                alpha(glass1, .1); 
                glass2 = ctmr_gauss_plot_edited(Lcrtx, [], [], S.cax, 0, S.cm, S.gsp); 
                alpha(glass2, .1); 
                plot3(em(depthch,1),em(depthch,2),em(depthch,3),'k.','markersize',10-5*(1/nch*10))
                if ~pltshowplanes(j)
                    plot3(em(nns,1), em(nns,2), em(nns,3), 'k.', 'markersize', 10-5*(1/nch*10)); 
                end 
            else 
                plot3(em(nns,1),em(nns,2),em(nns,3),'k.','markersize',10-5*(1/nch*10)) %plot electrodes
            end
            cameratoolbar('setmode',''); 
            litebrain(viewangle{j},.9); 
            wb=strcmpi(srf,'wholebrain'); 
            if any(wb)
                alpha(glass1,srfalpha{wb}); 
                alpha(glass2,srfalpha{wb}); 
            end
            if S.sliceplane=='c'
                view(180,270);
            elseif S.sliceplane=='a'
                view(180,0);
            end
            axis(axislim); 
            zoom(pltzoom(j)); 
            hold on; 
            colormap(gca,S.cm); 
            set(gca,'Clipping','off');
            % clone the surface plot for use in grab_slice
            surffig = figure(2);
            % set(gcf, 'visible', 'off');
            newhandle = copyobj(h, surffig);
            zoom(3.6);
            newhandle.Position(1:2) = [.28 .45];
            figure(1);
            clear srfplot
          case 'DEPTH' %plot depth electrode with parallel slice (plus surface behind it)
            eN=depths{j}; 
            [eNID,~,~]=intersect(find(nns),eN); %Get the specific channels for this depth, ignoring bad channels
            if isempty(eNID)
                axis off; 
                drows(drows==j)=[]; 
            elseif ~isempty(eNID)
                I.em=em; 
                I.w8s=w8s; 
                I.nns=nns; 
                sliceinfo(j).depthlabels=depthlabels{j};
                if S.sliceplane=='c'
                    OPSCEAsurfslice(pt,S.sliceplane,em(eNID,:),zeros(size(em(eNID,:))),opsceadatapath,[],S.cax,S.cm,S.gsp,j,1,force_angle_coronal(j));
                    if strcmp(sliceinfo(j).final_orientation, 'oc')
                        figure(2);
                        view(270,0); % flip surface brain to sagittal 
                        figure(1);
                    end
                elseif S.sliceplane=='a'
                    OPSCEAsurfslice_axial(pt,S.sliceplane,em(eNID,:),zeros(size(em(eNID,:))),opsceadatapath,[],S.cax,S.cm,S.gsp,j,1,force_angle_axial(j));
                end
                cameratoolbar('setmode','')
                axis off; 
                axis equal; 
                hold off; 
          
                %add light sources and set camera view
                delete(findall(gca,'type','light')); 
                view(sliceinfo(j).azel); %head-on angle for camlight reference
            
                hold on; % add color-coded planes parallel to
                         % slice, to highlight the plane of view
                         % during initial rotation below
                if S.sliceplane=='c'
                    camlight(0,30) %camlight at head-on angle and 25 degrees above azimuth (otherwise light reflects and will obscure cam view)
                    fill3(sliceinfo(j).corners(1,:), sliceinfo(j).corners(2,:)+.1,sliceinfo(j).corners(3,:),depthcolor{j},'edgecolor',depthcolor{j},'facealpha',0,'edgealpha',.5,'linewidth',3); 
                elseif S.sliceplane=='a'
                    % roll the brain
                    if sliceinfo(j).viewangle > 0 && sliceinfo(j).viewangle <= 90
                        % 1st quad, +
                        ang = -sliceinfo(j).viewangle;
                    elseif sliceinfo(j).viewangle > 90 && sliceinfo(j).viewangle <= 180
                        % 2nd quad, +
                        ang = 180 - sliceinfo(j).viewangle;
                    elseif sliceinfo(j).viewangle > 180 && sliceinfo(j).viewangle <= 270
                        % 3rd quad, +
                        ang = 180 - sliceinfo(j).viewangle;
                    elseif sliceinfo(j).viewangle > 270 && sliceinfo(j).viewangle <= 360
                        % 4th quad, +
                        ang = 360 - sliceinfo(j).viewangle;
                    elseif sliceinfo(j).viewangle > -90 && sliceinfo(j).viewangle <= 0
                        % 4th quad, -
                        ang = -sliceinfo(j).viewangle;
                    elseif sliceinfo(j).viewangle > -180 && sliceinfo(j).viewangle <= -90
                        % 3rd quad, -
                        ang = -(sliceinfo(j).viewangle + 180);
                    elseif sliceinfo(j).viewangle > -270 && sliceinfo(j).viewangle <= -180
                        % 2nd quad, -
                        ang = -(sliceinfo(j).viewangle + 180);
                    elseif sliceinfo(j).viewangle > -360 && sliceinfo(j).viewangle <= -270
                        % 1st quad, -
                        ang = -(sliceinfo(j).viewangle + 360);                        
                    end
                    
                    if sliceinfo(j).sagittal
                        if sliceinfo(j).viewangle < 90
                            ang = ang + 180;
                        end
                        camorbit(ang, 0, 'coordsys', [0 1 0]);
                        camroll(90);
                    else
                        camorbit(ang, 0, 'coordsys', [0 1 0]);
                    end
                    
                    if sliceinfo(j).sagittal
                        ang1 = 45;
                    else
                        angl = sliceinfo(j).viewangle;                    
                        if angl > 0 && angl <= 90
                            ang1 = -angl;
                        elseif angl > 90 && angl <= 180
                            ang1 = 180 - angl;
                            if ang1 < 15
                                ang1 = -ang1;
                            end
                        else
                            ang1 = 0;
                        end
                    end                    
                    
                    if sliceinfo(j).b < -40
                        ang2 = 30;
                    elseif sliceinfo(j).b >= -40 && sliceinfo(j).b < -15
                        ang2 = 40;
                    else
                        ang2 = 45;
                    end

                    camlight(ang1, ang2);
                    fill3(sliceinfo(j).corners(1,:), sliceinfo(j).corners(2,:),sliceinfo(j).corners(3,:)-.1,depthcolor{j},'edgecolor',depthcolor{j},'facealpha',0,'edgealpha',.5,'linewidth',3); 
                end                    
                hold off; 
                
                for showp=find(pltshowplanes) % now add color-coded slices planes on any subplot(s) where you indicated showplanes=1
                    if ~isempty(showp)
                        subplot(subplotrow(showp),subplotcolumn(showp),subplotnum{showp}); 
                        fill3(sliceinfo(j).corners(1,:), sliceinfo(j).corners(2,:), sliceinfo(j).corners(3,:), depthcolor{j},'edgecolor',depthcolor{j},'facealpha',.1,'edgealpha',.5,'linewidth',3);
                        subplot(subplotrow(j),subplotcolumn(j),subplotnum{j}); %switch back to the slice subplot at hand
                    end
                end
                axis(axislim); 
                camzoom(pltzoom(j)); % apply the specified zoom for that this view (usually similar for all depths but can depend on angle of slice, position of electrodes, etc)
                                
                if showlabels
                    ttl=depthlabels{j}; 
                else 
                    ttl=['Depth ' num2str(j-(find(isdepth==1,1)-1))]; 
                end

                ttloffset = -1; % vertical offset of title from bottom of axis, in millimeters
                if S.sliceplane=='c'
                    dim1 = mean(sliceinfo(j).corners(1,:));
                    dim2 = mean(sliceinfo(j).corners(2,:));
                    dim3 = axislim(5) + ttloffset;
                elseif S.sliceplane=='a'
                    if sliceinfo(j).sagittal                   
                        dim1 = mean(sliceinfo(j).corners(1,:));
                        dim2 = mean(sliceinfo(j).corners(2,:));
                        dim3 = axislim(5) + ttloffset;
                    else
                        dim1 = mean(sliceinfo(j).corners(1,:));
                        dim2 = axislim(3) + ttloffset;
                        dim3 = mean(sliceinfo(j).corners(3,:));
                    end
                end
                text(dim1, dim2, dim3, ttl, 'HorizontalAlignment', 'center', ...
                     'VerticalAlignment', 'top', 'color', depthcolor{j}, ...
                     'fontweight', 'bold', 'fontsize', 14);

                colormap(gca,S.cm);
                grab_slice(h, ptpath, ttl, j, depthcolor{j}, S.sliceplane, em, eNID);
            end
        end
    end
    close(2);
end

if dopdf
    close all;
    if ~exist('selected_leads', 'var')
        export_bs_figs(pt, channel_mat, [ptpath 'Imaging/Recon/figs']);
    end
    system(['python make_slice_pdf.py ' pt]);
end

function grab_slice(h, ptpath, label, n, color, orientation, em, eNID)
global sliceinfo;
targetdir = [ptpath 'Imaging/Recon/figs']; 
if ~exist(targetdir, 'dir')
    mkdir(targetdir);
end
newfig = figure; %('visible', 'off');
newhandle = copyobj(h, newfig); 
newhandle.Position(1:2) = [.28 .45];
if orientation=='c'
    zoomval = 4;
elseif orientation=='a'
    zoomval = 2;
end
camzoom(zoomval); % need to use camzoom instead of zoom after camorbit

fname = [sprintf('%02d', n) '_' label '_' orientation];

% managing dot size for electrodes
% 1: why does size of axial cuts fluctuate?
% 2: there must be a better way to get the image dimensions
% first write
exportgraphics(newfig, [targetdir '/' fname '.png']);
im = imread([targetdir '/' fname '.png']);
xdim = size(im, 1);
scalefactor = xdim/560; % most coronal cuts have width 560

figure(newfig);
camzoom(1/scalefactor);
hold on;
if orientation=='c'
    dotsize = round(14 - 6*(1-scalefactor));
    if sliceinfo(n).sagittal
        yoffset = 3;
    else
        yoffset = 1;
    end
    plot3(em(eNID,1), em(eNID,2)+yoffset, em(eNID,3), 'k-'); % depth probe (line between electrodes)
    plot3(em(eNID,1), em(eNID,2)+yoffset, em(eNID,3), 'k.','markersize',dotsize); % depth electrodes (dots)
elseif orientation=='a'
    dotsize = round(12 - 8*(1-scalefactor));
    if sliceinfo(n).sagittal
        plot3(em(eNID,1)-1, em(eNID,2), em(eNID,3), 'k-');
        plot3(em(eNID,1)-1, em(eNID,2), em(eNID,3), 'k.','markersize',dotsize); 
    else
    plot3(em(eNID,1), em(eNID,2), em(eNID,3)-1, 'k-');
    plot3(em(eNID,1), em(eNID,2), em(eNID,3)-1, 'k.','markersize',dotsize); 
    end
end
hold off;

% 2nd write
exportgraphics(newfig, [targetdir '/' fname '.png']);
close(newfig);

figure(2);
% set(gcf, 'visible', 'off');
fh = fill3(sliceinfo(n).corners(1,:), sliceinfo(n).corners(2,:), sliceinfo(n).corners(3,:), color,'edgecolor',color,'facealpha',.1,'edgealpha',.5,'linewidth',3);
fname = [sprintf('%02d', n) '_' label '_surf_' orientation];
exportgraphics(gcf, [targetdir '/' fname '.png']);
delete(fh);
if strcmp(sliceinfo(n).final_orientation, 'oc')
    view(180,270); % reset to default view for coronal
end
figure(1);

function ml = get_matching_labels(tablecol, arr)
ml = 1:3;
for i=1:length(arr)
    for j=1:length(tablecol)
        if strcmp(arr{i}, tablecol{j})
            ml = [ml j];
            break;
        end
    end
end
