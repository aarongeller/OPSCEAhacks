function newpial = undo_fs_ahectomy(pial, hipp, amyg)
% Glues hippocampus and amygdala back onto the MTL after they are
% removed by freesurfer.
%
% Requires one of:
% Automated Driving Toolbox  
% Sensor Fusion and Tracking Toolbox 
% UAV Toolbox  

pmesh = extendedObjectMesh(pial.vert, cast(pial.tri, 'double'));
hmesh = extendedObjectMesh(hipp.vert, cast(hipp.tri, 'double'));
amesh = extendedObjectMesh(amyg.vert, cast(amyg.tri, 'double'));

hamesh = join(hmesh, amesh);
phamesh = join(pmesh, hamesh);

newpial = [];
newpial.tri = phamesh.Faces;
newpial.vert = phamesh.Vertices;

