function [m,b,theta,e1]=get_mb(elecs, dimvec)
if size(elecs,1)>6
    e1=2; 
    e2=size(elecs,1)-1; 
else 
    e1=1; 
    e2=size(elecs,1); 
end % get the 2nd-the-last from each side (unless 6 or less contacts), helps in case offset or bent shaft
dims_chosen = find(dimvec);
dim1 = dims_chosen(1);
dim2 = dims_chosen(2);
m = (elecs(e1,dim2) - elecs(e2,dim2))/(elecs(e1,dim1) - elecs(e2,dim1));              
b = elecs(e1,dim2) - m*elecs(e1,dim1); %algebra: b = y-mx
centeringvec = zeros(1,3);
centeringvec(dim2) = b;
centered_elecs = elecs - centeringvec;
[thetas, rhos] = cart2pol(centered_elecs(:,dim1), centered_elecs(:,dim2));
theta = circ_mean(thetas(e1:e2));
m = tan(theta);
