function box=boxconstruct_2D(form,varargin)
% Construct a uniform box structure for different kinds of boxes, including
% Circle, AABB, OBB, POLYGON. Deal with only one box every time, so the
% input cannot be a cell array.

% Input:
%     form: 'Circle', 'AABB', 'OBB', 'POLYGON'.
%     varargin: Different information corresponding to differet kinds of boxes.
% Output:
%     box: The uniform structure of a box. For a Circle box,
%       box.center/radius. For other boxes, box.points/sa/id. 
%       The first and last element are NOT the same in box.points, but the boundary-character
%       is represented by box.id. And box.sa is a matrix,where the number
%       of ROWS equal number of separated axes.

% 2 kinds of id:
% In this function, box.points are all the boundary vertices and box.id is
% [1,2;2,3;...]; In function nrbcrvbox_2D.m, cpolyg are the BOUNDARY
% vertices of each box (comparing with the CONTROL vertices coefs of each box),
% and idtemcoefs means indices of cpolyg in coefs, for each box.

% box.points obtained by nrbcrvbox_2D.m: The first element is the same as
% the last one.
% box.points obtained by boxconstruct_2D.m: The first element is NOT the same as
% the last one.


box.form=form;

switch form
    case 'Circle' % varargin{1}=center, varargin{2}=r
        box.center=varargin{1};
        box.radius=varargin{2};
    case 'AABB' % varargin{1}=pt1=[Xmin,Ymin], varargin{2}=pt2=[Xmax,Ymax]
        pt1=varargin{1};pt1=pt1(:)';
        pt2=varargin{2};pt2=pt2(:)';
        box.points=[pt1;
            pt2(1),pt1(y);
            pt2;
            pt1(1),pt2(2)]; % size=[4,2]
        box.id=[1,2;2,3;3,4;4,1]; % The boundary-character is represented by the index instead of repeated points.
        box.sa=[1,0;0,1];
    case 'OBB'  % varargin{1}=vert, size=[2,4], varargin{2}=sa
        vert=varargin{1};
        sa=varargin{2};
        box.points=vert';
        box.id=[1,2;2,3;3,4;4,1];% In vert, the 4 points have been sorted.
        box.sa=sa';
    case 'POLYGON' % varargin{1}=cpolyg{i}, one polygon that may be obtained by nrbcrvbox_2D.
        cpolyg=varargin{1};
        idtemcoefs=varargin{2};
        num=size(cpolyg,1);
        if num==2 % AABB style, 'Interval'
            pt1=cpolyg(:,1);pt1=pt1';
            pt2=cpolyg(:,2);pt2=pt2';            
            box.points=[pt1;
            pt2(1),pt1(y);
            pt2;
            pt1(1),pt2(2)]; % size=[4,2]
            box.id=[1,2;2,3;3,4;4,1];
            box.sa=[1,0;0,1];
        else % NURBS control polygons, 'Control'
            % Calculate the separate axes of each edge of the polygon,which is the unit normal vector of each edge.            
            cpolyg_=cpolyg(2:end,:);
            cpolyg(end,:)=[];% size=[numpoints,2], and in cpolyg, the first element is the same as the last one
            edgevec=cpolyg_-cpolyg;
            edgevec=edgevec';
            edgevec=vecnorm(edgevec);
            normvec=[-edgevec(2,:);
            edgevec(1,:)];
            box.sa=normvec';
            box.points=cpolyg;
            box.id=[idtemcoefs(1:end-1),idtemcoefs(2:end)];
        end
    otherwise
        warning('The input box is invalid');
end
            
        
        
        
        
        
        
        
        
        
        
        