function [pt,sortuv]=pointinsert(pt1,pt2,pt3)
% Insert the 3rd point into a line connected by points pt1 and pt2 using
% the symbol of dot product.The direction of the line is from pt1 to pt2.

% Input:
%   pt1,pt2: 2D/3D coordinates of two points which are connected as a line.
%   pt3: The 3rd point which needs to be inserted into line pt1-pt2.

% Output:
%   pt: Coordinates of the total three points after inserting pt3, which
%       may be [pt1;pt2;pt3]/[pt1;pt3;pt2]/[pt3;pt1;pt2].
%   sortuv: Order of indices of pt1,pt2,pt3.

vector12=pt2-pt1;
vector13=pt3-pt1;
dot12_13=dot(vector12,vector13);
if (dot12_13<=0)
    pt=[pt3;pt1;pt2];
    sortuv=[3,1,2];
else
    vector23=pt3-pt2;
    dot21_23=dot(-vector12,vector23);
    if (dot21_23<=0)
        pt=[pt1;pt2;pt3];
        sortuv=[1,2,3];
    else
        pt=[pt1;pt3;pt2];
        sortuv=[1,3,2];
    end
end

end


