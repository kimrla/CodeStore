function box=boxconstruct_3D(form,varargin)

% Construct a uniform box structure for different kinds of boxes, including
% Sphere, AABB, OBB, POLYHEDRON. Deal with only one box every time, so the
% input cannot be a cell array.

% Input:
%     form: 'Sphere', 'AABB', 'OBB', 'POLYHEDRON'.
%     varargin: Different information corresponding to differet kinds of 
%       boxes, which are obtained by CircleBox.m, AABB_3D.m, OBB_3D.m,
%       nrbboundary.m.
% Output:
%     box: The uniform structure of a box. For a Sphere, box.center/radius. 
%       For other boxes, box.points/sa/idedge/idface/form. 
%       idedge stores index of points in each edge, and idface face stores
%       index of points in each rectangular or triangular face. 
%       box.sa is a matrix,where the number of ROWS equal number of all the
%       separated axes.
%       In box.points, all the points are vertices of the boundary box,
%       comparing with the control points in boundbox.points in  nrbboundary.m

% Notice: 1.Whether the 1st element is the same as the last one in
% box.points and box.id; 2. Whether the number of ROWS or COLUMNS equal the
% number of separated axes in box.sa.

% For Sphere:
% varargin{1}=center, varargin{2}=r
% For AABB:
% varargin{1}=pt1=[Xmin,Ymin, Zmin], varargin{2}=pt2=[Xmax,Ymax, Zmax]
% For OBB:
% varargin{1}=vert, size(vert)=[8,3], varargin{2}=sa
% For POLYHEDRON:
% varargin{1}=boundbox{i},obtained by nrbboundary.m
box.form=form;
switch form
    case 'Sphere' % varargin{1}=center, varargin{2}=r
        box.center=varargin{1};
        box.center=box.center(:)';% Row vector.
        box.radius=varargin{2};
    case 'AABB' % varargin{1}=pt1=[Xmin,Ymin, Zmin], varargin{2}=pt2=[Xmax,Ymax, Zmax]
        pt1=varargin{1};pt1=pt1(:)';
        pt2=varargin{2};pt2=pt2(:)';
        % 8 vertices of the box
        box.points=[pt1(1) pt1(2) pt1(3)
                    pt2(1) pt1(2) pt1(3)
                    pt2(1) pt2(2) pt1(3)
                    pt1(1) pt2(2) pt1(3)
                    pt1(1) pt1(2) pt2(3)
                    pt2(1) pt1(2) pt2(3)
                    pt2(1) pt2(2) pt2(3)
                    pt1(1) pt2(2) pt2(3)];
        box.idedge=[1,2
                    2,3
                    3,4
                    4,1
                    5,6
                    6,7
                    7,8
                    8,5
                    1,5
                    2,6
                    3,7
                    4,8];
        box.idface=[1,2,3,4;
                    5,6,7,8;
                    1,2,6,5;
                    2,3,7,6;
                    3,4,8,7;
                    4,1,5,8;];           
        box.sa=[1,0,0;0,1,0;0,0,1];% Number of ROWS equal number of sa       
    case 'OBB'  
        % varargin{1}=vert, size(vert)=[8,3], varargin{2}=sa(The number of
        % COLUMNS equal the number of sa).
        vert=varargin{1};
        sa=varargin{2};
        box.points=vert;
        box.idedge=[1,2
                    2,3
                    3,4
                    4,1
                    5,6
                    6,7
                    7,8
                    8,5
                    1,5
                    2,6
                    3,7
                    4,8];   
        box.idface=[1,2,3,4;
                    5,6,7,8;
                    1,2,6,5;
                    2,3,7,6;
                    3,4,8,7;
                    4,1,5,8;];  
        box.sa=sa';% The number of ROWS equal the number of sa
    case 'POLYHEDRON' % varargin{1}=boundbox{i},obtained by nrbboundary.m
        % Only process one box every time, meaning the input is not a cell
        % arry, but a structure array as boundbox{i}.
        boundbox=varargin{1};
        points=boundbox.points;
        idedge=[boundbox.id(:,1),boundbox.id(:,2)];
        idedge=[idedge;[boundbox.id(:,3),boundbox.id(:,2)]];
        idedge=[idedge;[boundbox.id(:,3),boundbox.id(:,1)]];
        idedge=sort(idedge,2);
        box.idedge=unique(idedge,'rows');        
        id=unique(boundbox.id);        
        box.points=points(id,:);
        h=@(x)(find(id==x));
        boundbox.idface=arrayfun(h,boundbox.id);      
        box.sa=boundbox.sa;       
    otherwise
        warning('The input box is invalid');
     
end
            