% DM = DistanceMatrix(datasites,centers)
%计算s维空间R中两个矩阵的距离
% Forms the distance matrix of two sets of points in R^s,
% i.e., DM(i,j) = || datasite_i - center_j ||_2.
%
% Input 
%   datasites: Nxs matrix representing a set of N data sites in R^s 
%              (i.e., each row contains one s-dimensional point)
%
%   centers:   Mxs matrix representing a set of M centers for RBFs in R^s
%              (also one center per row)
%
% Output
%   DM:        NxM matrix that contains the Euclidean distance between the i-th data site 
%              and the j-th center in the i,j position
%
% no error checking

function DM = DistanceMatrix(datasites,centers)

[N,s] = size(datasites);
[M,s] = size(centers);

% Accumulate sum of squares of coordinate differences
% The meshgrid command produces two NxM matrices:
%   center_cols consisting of N identical rows (containing the M centers)
%   data_rows consisting of M identical columns (containing the N data sites)
DM = zeros(N,M);

for d=1:s
    [center_cols,data_rows] = meshgrid(centers(:,d),datasites(:,d));
    DM = DM + (center_cols-data_rows).^2;
end
DM = sqrt(DM);



