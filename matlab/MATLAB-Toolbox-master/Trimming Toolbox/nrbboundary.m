function boundbox=nrbboundary(model,varargin)

% Construct the boundary of a NURBS or Bezier model(curve or surface).The
% boundary box is a polygon(2D curve) or polyhedron(3D curve or surface).
% For 2D, the pyhsical coords are [x,y] ; For 3D, the coords are [x,y,z].

% Input:
%   model: NURBS structure nrb or Bezier structure bzr.
%   varargin{1}='whole': Indicate that the input is NURBS model and the
%       output box is the whole poly-box of  the whole model without seg.
% Output:
%   boundbox: The boundary box constructed by the control points, which is
%       a structure(for Bezier model) or a structure-cell array (for NURBS model). 
%       Refer to boxconstruct_2D('POLYGON'), boundary.m and convhull.m.

% The box's structure.
%   boundbox.points: Physical coordinates (2D or 3D) of the vertices of all 
%       the control points (not must be the box's vertices).  If the box is 
%       2D, the first point is NOT the same as the last one and size(points)=[dim,2]. 
%       In each cell, the number of ROWS equals to the number of control points.
%   boundbox.id: Index of the box's vertices, where id=[dim,2] for 2D polygon,
%       storing the 2 vertices-id of each edge; and id=[dim,3] for 3D polyhedron, 
%       storing the 3 vertices-id of each triangular patch. 
%       If boundbox.points is a cell array, then boundbox.id is also a corresponding cell array. 
%   boundbox.sa: Separated axes, which is the unit normal vectors of each edge or triangular patch.
%       The number of ROWS equal the number of edges or patches in each box.
%       The sequence corresponds to boundbox.id. 

points=model.coefs;
dimdet=sum(points(3,:));
knots=model.knots;
ord=model.order;
nvec=@(p1,p2,p3)(vecnorm(cross((p2-p1)',(p3-p2)')));
% Calculate unit normal vector of a 3D triangular patch using the 3 vertices' coords.
% If cross(A,B), A,B are matrixes, then calculate cross-products of each column of A and B.
if strcmp(model.form,'B-NURBS')
    if (length(ord)==1 && dimdet<eps) % 2D NURBS curve.The box is polygons.
        % Construct default number of boxes equals to num(break-points)-1, by using knot insert.
        [cpolyg,idtem] = nrbcrvbox_2D( model,'Control'); 
        % Cell array. The forst point is the same as the last one.
        numbox=length(cpolyg);
        % construct box.id,box.sa and remove the last repeated point of each cpolygon.
        for i=1:numbox           
        % box.points obtained by nrbcrvbox_2D.m: The first element is the same as the last one.
        % box.points obtained by boxconstruct_2D.m: The first element is NOT the same as the last one.            
            boundbox{i}=boxconstruct_2D('POLYGON',cpolyg{i},idtem{i}); 
        end
  
    else % 3D NURBS curve or surface. The box is polyhedrons, enveloped by triangular patches.
        if nargin==1% Default: Transform the NURBS to Bezier
            bzr=nrb2bzr(model);% Transform the NURBS curve or surface to Bezier model.
            numbzr=numel(bzr);
            for i=1:numbzr
                pts=bzr{i}.coefs;
                pts=pts(1:3,:)./repmat(pts(4,:),3,1);pts=pts';
                id=convhull(pts(:,1),pts(:,2),pts(:,3));
                boundbox{i}.points=pts;
                boundbox{i}.id=id;
                % Calculate the unit normal vector of each triangular patch.
                temsa=nvec(pts(id(:,1),:),pts(id(:,2),:),pts(id(:,3),:));
                boundbox{i}.sa=temsa';% ROWS equal to the number of triangular patches.
            end         
        else
            if strcmp(varargin{1},'whole')
                pts=model.coefs;
                pts=pts(1:3,:)./repmat(pts(4,:),3,1);pts=pts';
                id=convhull(pts(:,1),pts(:,2),pts(:,3));
                boundbox.points=pts;
                boundbox.id=id;
                temsa=nvec(pts(id(:,1),:),pts(id(:,2),:),pts(id(:,3),:));
                boundbox.sa=temsa';
            end
        end
    end
          
elseif strcmp(model.form,'Bezier')
    pts=model.coefs;
    pts=pts(1:3,:)./repmat(pts(4,:),3,1);
    pts=pts';
    id=convhull(pts(:,1),pts(:,2),pts(:,3));
    boundbox.points=pts;
    boundbox.id=id;
    % Calculate the unit normal vector of each triangular patch.
    temsa=nvec(pts(id(:,1),:),pts(id(:,2),:),pts(id(:,3),:));
    boundbox.sa=temsa';% ROWS equal to the number of triangular patches.
end

end

%% demo
% srf=nrbtestsrf;
% box=nrbboundary(srf);
% figure;nrbctrlplot(srf)
% hold on;
% k=3;
% trimesh(box{k}.id,box{k}.points(:,1),box{k}.points(:,2),box{k}.points(:,3));
% k=7;
% trimesh(box{k}.id,box{k}.points(:,1),box{k}.points(:,2),box{k}.points(:,3));



