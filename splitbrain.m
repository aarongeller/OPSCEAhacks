function split = splitbrain(cortex, orientation, b, m)
% Splits the brain and exports a mesh that is a subset of the original
% mesh. Requires the function splitFV (Matlab Exchange)
%
%
% INPUTS:
% cortex - the original mesh. a structure containing tri and vert 
% 
% orientation - slice orientation ('c', 'a', 's', or 'oc'), with
% slices defined as follows, for axes defined medial/lateral = X
% axis, anterior/posterior = Y axis, inferior/superior = Z axis:
%
% c (coronal): defined by a line in the XY plane, angle near n*pi
%              for n in integers including zero
% a (axial): defined by a line in the XZ plane
% s (sagittal): defined by a line in the XZ plane, angle near
%               n*pi/2 for n in integers excluding zero
% oc (oblique coronal): defined by a line in the YZ plane
%
% b- intercept
%
% m- slope
%
% OUTPUTS:
% split - a new mesh split at the specified point. A structure containing
% tri and vert

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


%mind the gap (in mm): prevents triangles that span the gap from inducing unwanted mesh
thegap = max([abs(5.*m) 10]);

if strcmp(orientation, 'c')
    if (min(cortex.cortex.vert(:,2)) < b) && (b < max(cortex.cortex.vert(:,2)))
        % coronal view uses a cut in the XY plane
        [idx, ~] = sort(find(abs((cortex.cortex.vert(:,2) - (m.*cortex.cortex.vert(:,1) + b + thegap))) <= thegap)); 
        % get indices of verts with y between (mx + b) and (mx + b + 2*thegap)
    else 
        display('splitbrain.m: intercept lies outside A-P extent of the mesh.');
        [idx, ~] = sort(find(abs((cortex.cortex.vert(:,2) - (m.*cortex.cortex.vert(:,1) + b + thegap))+2*thegap)<=thegap));
    end
    
elseif strcmp(orientation, 'a') || strcmp(orientation, 's')
    % both axial and sagittal views use a cut in the XZ plane
    if (min(cortex.cortex.vert(:,3)) < b) && (b < max(cortex.cortex.vert(:,3)))
        [idx, ~] = sort(find(abs((cortex.cortex.vert(:,3) - (m.*cortex.cortex.vert(:,1) + b - thegap))) <= thegap)); 
        % get indices of verts with z between (mx + b - 2*thegap) and (mx + b)
    else 
        display('splitbrain.m: intercept lies outside S-I extent of the mesh.');
        [idx, ~] = sort(find(abs((cortex.cortex.vert(:,3)-(m.*cortex.cortex.vert(:,1) + b + thegap))+2*thegap)<=thegap));
    end

elseif strcmp(orientation, 'oc')
    % oblique coronal view uses a cut in the YZ plane
    [idx, ~] = sort(find(abs((cortex.cortex.vert(:,3) - (m.*cortex.cortex.vert(:,2) + b + thegap))) <= thegap)); 
    % get indices of verts with z between (my + b + 2*thegap) and (my + b)
end

split.vert = [];
split.tri = [];

if isempty(idx) && ~strcmp(orientation, 's')
    % if we're in sagittal mode and there is no split, then we are working
    % with the L hemisphere and need to *omit* the whole hemisphere;
    % otherwise include all of it
    split.vert = cortex.cortex.vert;
    split.tri = cortex.cortex.tri;
elseif ~isempty(idx)
    mesh.tri = delete_verts(cortex.cortex.tri, idx);        
    mesh.vert = cortex.cortex.vert;

    disp('Generating partial mesh for slice view...')
    FVout = splitFV(mesh.tri,mesh.vert);

    sv = checkFV(FVout);
    [~,svinds] = sort(sv);
    foundgood = 0;
    for i=length(svinds):-1:1
        maxmeshind = svinds(i);
        fv.vert = FVout(maxmeshind).vertices;
        fv.tri = FVout(maxmeshind).faces;

        if orientation_good(fv.vert, m, b, orientation)
            foundgood = 1;
            break;
        elseif i==1
            % no submesh looked good
            break;
        end
        display(['splitFV chose wrong submesh, trying submesh #' int2str(svinds(i-1)) '...']);
    end

    if foundgood
        split = fv;
    elseif ~isempty(idx)
        beep; pause(0.25); beep;
        meshind = input('*** Automated submesh selection failed! Enter submesh to use. *** ');
        split.vert = FVout(meshind).vertices;
        split.tri = FVout(meshind).faces;
    end
end

function tri = delete_verts(tri, idx)
while(~isempty(intersect(tri(:,2), idx)))
    [~,ia2,~] = intersect(tri(:,2), idx); %looks for the tris that match the indices of the verts on the line
    tri(ia2,:) = []; %remove the tris from the mesh dataset that have already been counted
end

function status = orientation_good(verts, m, b, orientation)
% check that the chunk we're choosing is on the right (correct) side
% of the line.

status = 0;
centroid = mean(verts);

if strcmp(orientation, 'c') && centroid(2) < m*centroid(1) + b
    % line is in XY plane
    % for coronal cut want to choose the part behind the plane (we're
    % looking back) so we want centroid below the line
    status = 1;
elseif strcmp(orientation, 'a') && centroid(3) > m*centroid(1) + b
    % line is in XZ plane
    % for axial cut we want to choose the part above the plane
    % (we're looking up) so we want centroid above the line
    status = 1;
elseif strcmp(orientation, 's') && centroid(1) > (centroid(3) - b)/m
    % line is in XZ plane
    % for sagittal cut we want the part to the right of the plane
    % (we're looking to the right, which has increasing X)
    status = 1;
elseif strcmp(orientation, 'oc') && centroid(2) < (centroid(3) - b)/m
    % line is in YZ plane
    % for oblique coronal cut we want the part to the left of the
    % plane when viewed from the side (i.e. the posterior part);
    % note sign flipped from conventional sagittal view
    status = 1;
end

function sizevec=checkFV(fvout)
sizevec = [];
for i=1:length(fvout)
    sizevec(end+1) = size(fvout(i).vertices,1);
end
