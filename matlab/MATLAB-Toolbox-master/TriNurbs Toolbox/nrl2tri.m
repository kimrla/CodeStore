function tnrl=nrl2tri(nurl, h0)

% nrb2tri: Transform a nurl geometry into triangular representation.
% 
% Calling Sequences:
% 
%     tnrl=nrl2tri(nurl, h0)
% 
% INPUTS:
%
%      nurl - A nurl curve, surface or volume.
%
%      h0 - Mesh seed length.
%
% OUTPUT:
% 
%     tnrl - Triangular representation of the nurl (tri-nurbs) surface.
%               A structure array with the following fields:
%
%         tnrl.form =  'Tri-NURL'.  Triangular representation 
%                               of nurbs surface
%
%         tnrl.surface - The nurl surface.
%
%         tnrl.numbers - The number of nodes or points
%                                  and the number of triangles.
%
%         tnrb.seeds    -  Mesh seed length of the surface and the parametric domain.
%
%         tnrl.nodes  -  Nodes on parametric domain.
%
%         tnrl.points  -  Points on the nurbs surface.
%
%         tnrl.delaunay  -  Delaunay triangular reprezetation.
%  

if (~isstruct(nurl))
  error('NURL representation is not structure!');
end

if (~strcmp(nurl.form,'L-NURL'))
  error('Not a recognised NURL representation');
end

if (iscell(nurl.knots))
    if (size(nurl.knots,2) == 3)
    %% NURL structure represents a volume
        
    elseif (size(nurl.knots,2) == 2)
    %% NURL structure represents a surface
        % Transpose the NURL surface to make the parametric curves 
        %   on u direction approximately with equal length.
        [Lu, Lv]=nrlsrfmeasure(nurl);
        ru=max(Lu)/min(Lu);
        rv=max(Lv)/min(Lv);
        if nargin==1
            h0=max([Lu, Lv])/30;
        end
        if ru<rv
            nurl=nrltransp(nurl);
            Lv=Lu;
        end

        % Delaunay Triangulation
        m=round(max(Lv)/h0); 
        nodes=cell(1, m); 
        s=linspace(0, 1, m);
        [Lu, ~]=nrlsrfmeasure(nurl, s);
        nd=0; num=zeros(m, 1);
        for i=1:m
            num(i)=max([ceil(Lu(i)/h0), 1]);
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
        pnts=nrleval(nurl, [x'; y']);

        % Make triangular representation of nurbs surface
        tnrl.form='Tri-NURL';
        tnrl.dim=2;
        tnrl.nurl=nurl;
        tnrl.numbers=[nd, nt];
        tnrl.seeds=[h0, 1/m];
        tnrl.nodes=[x, y];
        tnrl.points=pnts';
        tnrl.delaunay=sort(tri, 2);
    end
else
%% NURL structure represents a curve
    % Transform a nurl curve into line segment representation
    crvs=nrlsplits(nurl, nurl.intervals);
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
        nodes{i}=linspace(nurl.intervals(i), nurl.intervals(i+1), num(i));    
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
    pnts=nrleval(nurl, u);

    % Make triangular representation of nurbs surface
    tnrl.form='Tri-NURL';
    tnrl.dim=1;
    tnrl.nurbs=nurl;
    tnrl.numbers=[nd, nd-1];
    tnrl.seeds=[h0, 1/nd];
    tnrl.nodes=u';
    tnrl.points=pnts';
    tnrl.delaunay=[1:nd-1; 2:nd]';
end


%% demo - curve and plane
% % The mesh seed length (h0)
% h0=0.5;
% 
% % Create a plane square and a plane cuve
% lin1=nrlline([0,1], [9,1]);
% lin2=nrlline([0,6], [9,6]);
% srf=nrlruled(lin1, lin2);
% crv=nrltestcrv;
% 
% % Transform a nurbs surface into triangular representation
% tsrf=nrl2tri(srf, h0);
% tcrv=nrl2tri(crv, h0);
% 
% % Plot results
% figure; hold on; 
% tnrlplot(tsrf);
% tnrlplot(tcrv);
% axis equal;


%% demo - surface
% % The mesh seed length (h0)
% h0=0.2;
% 
% % Create a nurbs sphere
% circ=nrlcirc(1, [0,0], 0, pi);
% srf=nrlrevolve(circ, [0,0,0], [1,0,0], 2*pi);
% % srf=nrltestsrf;
% 
% % Transform a nurbs surface into triangular representation
% tnrl=nrl2tri(srf, h0);
% 
% % Plot the nodes on parametric domain
% figure; hold on;
% plot(tnrl.nodes(:,1), tnrl.nodes(:,2), 'ro');
% axis equal;
% title('Nodes on parametric domain');
% 
% % Plot the surface
% figure; hold on;
% nrlplot(srf, [35, 35]);
% view(3);
% shading interp;
% 
% figure; triplot(tnrl.delaunay, tnrl.nodes(:,1), tnrl.nodes(:,2)); 
% title('Parametric mesh');
% 
% figure; 
% trisurf(tnrl.delaunay, tnrl.points(:,1), tnrl.points(:,2), tnrl.points(:,3));
% axis equal;
% title('Geometric grid');







