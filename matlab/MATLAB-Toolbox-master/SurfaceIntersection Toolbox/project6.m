%% ========================================================================
function pass = project6( p )
% project all 6 vertices of both triangles onto vector p and check if two
% projections overlap
global V1 V2 V3 U1 U2 U3
dot_prod = @(a,b) a(:,1).*b(:,1)+a(:,2).*b(:,2)+a(:,3).*b(:,3);
%% Project vertices of triangle 1 and find the bounds min1 and max1
P = [dot_prod(p, V1), dot_prod(p, V2), dot_prod(p, V3)];
max1 = max(P,[],2);
min1 = min(P,[],2);
%% Project vertices of triangle 2 and find the bounds min1 and max1
P = [dot_prod(p, U1), dot_prod(p, U2), dot_prod(p, U3)];
max2 = max(P,[],2);
min2 = min(P,[],2);
%% Compare the bounds to see if they overlap
pass = (( min1 < max2 ) & ( min2 < max1 )) | ~dot_prod(p, p);
end