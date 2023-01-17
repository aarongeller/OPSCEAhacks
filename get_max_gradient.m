function dim=get_max_gradient(elecs)
dim = nan;
dx = abs(elecs(end, 1) - elecs(1, 1));
dy = abs(elecs(end, 2) - elecs(1, 2));
dz = abs(elecs(end, 3) - elecs(1, 3));

if dx > dy && dx > dz
    dim = 1;
elseif dy > dx && dy > dz
    dim = 2;
else
    dim = 3;
end
