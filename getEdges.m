function [xedge, yedge, zedge] = getEdges(alphamap, XX, YY, ZZ, orientation, theta)
%find the first row that has a nonzero
%
% Omni-planar and surface casting of epileptiform activity (OPSCEA)
% 
% Dr. Jon Kleen, 2017

if ~exist('orientation', 'var')
    orientation = 'c';
end

if ~exist('theta', 'var')
    theta = 0;
end

% if strcmp(orientation, 'oc')
%     if theta>=0
%         alphamap = imrotate(alphamap, -90);
%     else
%         alphamap = imrotate(alphamap, -270);
%     end
% end

i=1;
allblack = 1;
while allblack==1

    if size(find(alphamap(i,:)==1))>0
        firstrow = i;
        allblack = 0;
    end
    
    i = i+1;
end

%find the first column that has a nonzero
j=1;
allblack=1;
while allblack==1

    if size(find(alphamap(:,j)==1))>0
        firstcol = j;
        allblack = 0;
    end
    
    j = j+1;
end

%same thing coming from the bottom
i = size(alphamap,1);
allblack = 1;
while allblack==1

    if size(find(alphamap(i,:)==1))>0
        lastrow = i;
        allblack = 0;
    end
    
    i = i-1;
end

j=size(alphamap,2);
allblack=1;
while allblack==1

    if size(find(alphamap(:,j)==1))>0
        lastcol = j;
        allblack = 0;
    end
    
    j = j-1;
end

% corner1 = [firstrow, firstcol];
% corner2 = [firstrow, lastcol];
% corner3 = [lastrow, firstcol];
% corner4 = [lastrow, lastcol];

if strcmp(orientation, 'c')
    xedge = ([XX(firstrow, firstcol) XX(firstrow, lastcol)]);
    yedge = ([YY(firstrow, firstcol) YY(firstrow, lastcol)]);
    zedge = ([ZZ(lastrow, firstcol) ZZ(firstrow, firstcol)]);
elseif strcmp(orientation, 'a')
    xedge = ([XX(firstrow, firstcol) XX(firstrow, lastcol)]);
    yedge = ([YY(lastrow, firstcol) YY(firstrow, firstcol)]);
    zedge = ([ZZ(firstrow, firstcol) ZZ(firstrow, lastcol)]);    
elseif strcmp(orientation, 'oc')
    xedge = ([XX(lastrow, firstcol) XX(firstrow, firstcol)]); % like Z
    yedge = ([YY(firstrow, firstcol) YY(firstrow, lastcol)]); % like X
    zedge = ([ZZ(firstrow, firstcol) ZZ(firstrow, lastcol)]); % like Y
end
