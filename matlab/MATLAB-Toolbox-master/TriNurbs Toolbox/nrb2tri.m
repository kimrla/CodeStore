function tnrb=nrb2tri(nurbs, h0)

% nrb2tri: Transform a nurbs surface into triangular representation.
% 
% Calling Sequences:
% 
%     tnrb=nrb2tri(nurbs, h0)
% 
% INPUTS:
%
%      nurbs - A nurbs curve, surface or volume.
%
%      h0 - Mesh seed length.
%
% OUTPUT:
% 
%     tnrb - Triangular representation of the nurbs (tri-nurbs) surface.
%               A structure array with the following fields:
%
%         tnrb.form =  'Tri-NURBS'.  Triangular representation 
%                               of nurbs surface
%
%         tnrb.surface - The nurbs surface.
%
%         tnrb.numbers - The number of nodes or points
%                                  and the number of triangles.
%
%         tnrb.seeds    -  Mesh seed length of the surface and the parametric domain.
%
%         tnrb.nodes  -  Nodes on parametric domain.
%
%         tnrb.points  -  Points on the nurbs surface.
%
%         tnrb.delaunay  -  Delaunay triangular reprezetation.
%  

if (~isstruct(nurbs))
  error('NURBS representation is not structure!');
end

if (~strcmp(nurbs.form,'B-NURBS'))
  error('Not a recognised NURBS representation');
end

if (iscell(nurbs.knots))
    if (size(nurbs.knots,2) == 3)
    %% NURBS structure represents a volume
        
    elseif (size(nurbs.knots,2) == 2)
    %% NURBS structure represents a surface
        % Transpose the NURBS surface to make the parametric curves 
        %   on u direction approximately with equal length.
        [Lu, Lv]=nrbsrfmeasure(nurbs);
        ru=max(Lu)/min(Lu);
        rv=max(Lv)/min(Lv);
        if nargin==1
            h0=max([Lu, Lv])/30;
        end
        if ru<rv
            nurbs=nrbtransp(nurbs);
            Lv=Lu;
        end

        % Delaunay Triangulation
        m=round(max(Lv)/h0); 
        nodes=cell(1, m); 
        s=linspace(0, 1, m);
        [Lu, ~]=nrbsrfmeasure(nurbs, s);
        nd=0; num=zeros(m, 1);
        for i=1:m
            num(i)=max([ceil(Lu(i)/h0), 1]); % Essential!!!!!!
            nodes{i}=linspace(0, 1, num(i));
            nd=nd+num(i);
        end
        x=zeros(nd, 1); y=x; 
        t=0;
        for i=1:m
            ni=num(i); 
            x(t+1:t+ni)=s(i);
            y(t+1:t+ni)=nodes{i};
            t=t+ni;
        end
        t=0; nt=0;
        trgc=cell(1, m-1);
        nc=zeros(m-1, 1);
        for i=1:m-1
            ni=num(i)+num(i+1); 
            xi=x(t+1:t+ni);
            yi=y(t+1:t+ni);
            p=[xi, yi];
            trgi = delaunayTriangulation(p);
            trgc{i}=trgi(:,:)+t;
            t=t+num(i);    
            nc(i)=size(trgi, 1);
            nt=nt+nc(i);
        end
        tri=zeros(nt, 3);
        t=0;
        for i=1:m-1
            tri(t+1:t+nc(i), :)=trgc{i};
            t=t+nc(i);
        end
        pnts=nrbeval(nurbs, [x'; y']);

        % Make triangular representation of nurbs surface
        tnrb.form='Tri-NURBS';
        tnrb.dim=2;
        tnrb.nurbs=nurbs;
        tnrb.numbers=[nd, nt];
        tnrb.seeds=[h0, 1/m, 1/max(num)];
        tnrb.nodes=[x, y];
        tnrb.points=pnts';
        tnrb.delaunay=sort(tri, 2);
    end
else
%% NURBS structure represents a curve
    % Transform a nurbs curve into line segment representation
    nrl=nrb2nrl(nurbs);
    crvs=nrlsplits(nrl, nrl.intervals);
    n=length(crvs);
    L=zeros(1, n);
    for i=1:n
        L(i)=nrlmeasure(crvs(i));
    end
    nodes=cell(1, n); 
    nd=0; 
    num=zeros(n, 1);
    for i=1:n
        num(i)=max([ceil(L(i)/h0), 2]);
        nodes{i}=linspace(nrl.intervals(i), nrl.intervals(i+1), num(i));    
        nd=nd+num(i)-1;
    end
    nd=nd+1;
    u=zeros(1, nd);
    t=0;
    for i=1:n
        ni=num(i); 
        u(t+1:t+ni-1)=nodes{i}(1:ni-1);
        t=t+ni-1;
    end
    u(end)=nodes{i}(end);
    pnts=nrbeval(nurbs, u);

    % Make triangular representation of nurbs surface
    tnrb.form='Tri-NURBS';
    tnrb.dim=1;
    tnrb.nurbs=nurbs;
    tnrb.numbers=[nd, nd-1];
    tnrb.seeds=[h0, 1/nd];
    tnrb.nodes=u';
    tnrb.points=pnts';
    tnrb.delaunay=[1:nd-1; 2:nd]';
end


%% demo - curve and plane
% % The mesh seed length (h0)
% h0=0.5;
% 
% % Create a plane square and a plane cuve
% lin1=nrbline([0,1], [9,1]);
% lin2=nrbline([0,6], [9,6]);
% srf=nrbruled(lin1, lin2);
% crv=nrbtestcrv;
% 
% % Transform a nurbs surface into triangular representation
% tsrf=nrb2tri(srf, h0);
% tcrv=nrb2tri(crv, h0);
% 
% % Plot results
% figure; hold on; 
% tnrbplot(tsrf);
% tnrbplot(tcrv);
% axis equal;


%% demo - surface
% % The mesh seed length (h0)
% h0=0.2;
% 
% % Create a nurbs sphere
% circ=nrbcirc(1, [0,0], 0, pi);
% srf=nrbrevolve(circ, [0,0,0], [1,0,0], 2*pi);
% % srf=nrbtestsrf;
% 
% % Transform a nurbs surface into triangular representation
% tnrb=nrb2tri(srf, h0);
% 
% % Plot the nodes on parametric domain
% figure; hold on;
% plot(tnrb.nodes(:,1), tnrb.nodes(:,2), 'ro');
% axis equal;
% title('Nodes on parametric domain');
% 
% % Plot the surface
% figure; hold on;
% nrbplot(srf, [35, 35]);
% view(3);
% shading interp;
% 
% figure; triplot(tnrb.delaunay, tnrb.nodes(:,1), tnrb.nodes(:,2)); 
% title('Parametric mesh');
% 
% figure; 
% trisurf(tnrb.delaunay, tnrb.points(:,1), tnrb.points(:,2), tnrb.points(:,3));
% axis equal;
% title('Geometric grid');






