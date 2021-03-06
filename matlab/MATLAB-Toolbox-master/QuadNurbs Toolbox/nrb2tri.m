function tnrb=nrb2tri(nurbs, h0)

% nrb2tri: Transform a nurbs geometry into triangular representation.
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
%         tnrb.edges - A cell array of edge nodes in sequence.
%
%         tnrb.curves - A structure array of edge curves.
%
%         tnrb.cnodes - A cell array of parametric nodes of the edge curves.
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
        % Measure the length of parameric curves on a nurbs volume
        [Lu, Lv, Lw]=nrbvolmeasure(nurbs);
        L={Lu, Lv, Lw};
        rt(1)=max(Lu)/min(Lu);
        rt(2)=max(Lv)/min(Lv);
        rt(3)=max(Lw)/min(Lw);
        [~, order]=sort(rt);
        nurbs=nrbpermute(nurbs, order);

        % Get nodes on each layers of the volume
        Lu=L{order(1)};
        m=max([round(max(Lu)/h0)+1, 2]); 
        r=linspace(0, 1, m);
        tsrfs=cell(m, 1);
        nd=0; num=zeros(m, 1); 
        for i=1:m
            tsrfs{i}=nrl2tri(nrlvolextract(nurbs, r(i), [], []), h0);
            num(i)=tsrfs{i}.numbers(1);
            nd=nd+num(i);
        end
        pts=zeros(nd, 3); pnts=pts;
        t=0;
         for i=1:m
            ni=num(i); 
            pts(t+1:t+ni, 1)=r(i)*ones(num(i), 1);
            pts(t+1:t+ni, 2:3)=tsrfs{i}.nodes;
            pnts(t+1:t+ni, :)=tsrfs{i}.points;
            t=t+ni;
         end

        % Delaunay Triangulation
        t=0; nt=0;
        trgc=cell(1, m-1);
        nc=zeros(m-1, 1);
        for i=1:m-1
            ni=num(i)+num(i+1); 
            ptsi=pts(t+1:t+ni, :);
            trgi = delaunayTriangulation(ptsi);
            trgc{i}=trgi(:,:)+t;
            t=t+num(i);
            nc(i)=size(trgi, 1);
            nt=nt+nc(i);
        end
        tri=zeros(nt, 4);
        t=0;
        for i=1:m-1
            tri(t+1:t+nc(i), :)=trgc{i};
            t=t+nc(i);
        end

        % Make triangular representation of nurbs surface
        tnrb.form='Tri-NURBS';
        tnrb.dim=3;
        tnrb.nurbs=nurbs;
        tnrb.numbers=[nd, nt];
        tnrb.seeds=h0;
        tnrb.nodes=pts;
        tnrb.points=pnts;
        tnrb.delaunay=triangulation(tri, pnts);
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

        % Get nodes, curves and edge node numbers
        m=max([round(max(Lv)/h0)+1, 2]); 
        nodes=cell(1, m); 
        s=linspace(0, 1, m);
        [Lu, ~]=nrbsrfmeasure(nurbs, s);
        nd=0; num=zeros(m, 1);
        for i=1:m
            if Lu(i)<h0*1e-3
                num(i)=max([round(Lu(i)/h0), 1]);
            else
                num(i)=max([round(Lu(i)/h0)+1, 2]);
            end
            nodes{i}=linspace(0, 1, num(i));
            nd=nd+num(i);
        end
        x=zeros(nd, 1); y=x; 
        t=0; ed1=zeros(1,m); ed2=ed1;
        for i=1:m
            ni=num(i); 
            ed1(i)=t+1;
            x(t+1:t+ni)=s(i);
            y(t+1:t+ni)=nodes{i};
            t=t+ni;
            ed2(i)=t;
        end
        ed3=1:num(1);
        ed4=nd-num(end)+1:nd;
        crvs=nrbextract(nurbs);
        cnodes={linspace(0,1,m), linspace(0,1,m), ...
            linspace(0,1,num(1)), linspace(0,1,num(end))};
        
        % Delaunay Triangulation
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
        tnrb.seeds=h0;
        tnrb.nodes=[x, y];
        tnrb.points=pnts';
        tnrb.delaunay=tri;
        tnrb.edges={ed1, ed2, ed3, ed4}; 
        tnrb.curves=crvs([3,4,1,2]);
        tnrb.cnodes=cnodes;
        tnrb.pts2tri=tnrbpts2tri(tnrb);
        
        % Remove edges that are points
        id=true(1,4);
        for i=1:4
             if length(tnrb.edges{i})<=1
                 id(i)=false;
             end
        end
        tnrb.edges=tnrb.edges(id);
        tnrb.curves=tnrb.curves(id);
        tnrb.cnodes=tnrb.cnodes(id);
        tnrb.numbers(3)=length(tnrb.edges);
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
    tnrb.seeds=h0;
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


%% demo - surface 1
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


%% demo - surface 2
% % The mesh seed length (h0)
% h0=0.2;
% 
% % Create a nurbs sphere
% circ1=nrbcirc(1, [0,0], 0, pi);
% circ2=nrbcirc(1, [0,0], pi, 2*pi);
% circ2=nrbreverse(circ2);
% srf=nrbruled(circ1, circ2);
% figure; nrbplot(srf, [10,10]);
% view(2);
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
% figure; triplot(tnrb.delaunay, tnrb.nodes(:,1), tnrb.nodes(:,2)); 
% title('Parametric mesh');
% 
% figure; hold on;
% trisurf(tnrb.delaunay, tnrb.points(:,1), tnrb.points(:,2), tnrb.points(:,3));
% axis equal; view(2);
% title('Geometric grid');


%% demo - volume - 1
% % The mesh seed length (h0)
% h0=0.3;
% 
% % Creat a nurbs volume and transform to nurl 
% crv1 = nrbcirc(1,[0 0],0, pi/2);
% crv2 = nrbcirc(2,[0 0],0, pi/2);
% srf = nrbruled (crv1, crv2);
% srf = nrbtform (srf, [1 0 0 0; 0 1 0 1; 0 0 1 0; 0 0 0 1]);
% vol = nrbrevolve (srf, [0 0 0], [1 0 0], pi/2);
% figure; nrbplot(vol, [15, 15, 15], 'light', 'on');
% 
% % Transform the nurbs volume into triangular representation
% tnrb=nrb2tri(vol, h0);
% figure; tetramesh(tnrb.delaunay, tnrb.points);


%% demo - volume 2
% % The mesh seed length (h0)
% h0=0.3;
% 
% % Creat a nurbs volume and transform to nurl 
% % crv1=nrbline([0,0,0], [1,0,0]);
% % crv2=nrbline([0,1,0], [1,1,0]);
% % srf = nrbruled (crv1, crv2);
% % vol=nrbextrude(srf, [0,0,0.1]);
% % crv=nrbline([0,0,0], [1,0,0]);
% % srf=nrbrevolve(crv, [0,0,0], [0,0,1], pi/2);
% % vol=nrbextrude(srf, [0,0,0.1]);
% crv1 = nrbcirc(1,[0 0],0, pi/2);
% crv2 = nrbcirc(1.1,[0 0],0, pi/2);
% srf = nrbruled (crv1, crv2);
% srf = nrbtform (srf, [1 0 0 0; 0 1 0 1; 0 0 1 0; 0 0 0 1]);
% vol = nrbrevolve (srf, [0 0 0], [1 0 0], pi/2);
% figure; nrbplot(vol, [15, 15, 2]);
% 
% % Transform the nurbs volume into triangular representation
% tnrb=nrb2tri(vol, h0);
% points=tnrb.points;
% TR = triangulation(tnrb.delaunay, points);
% FBtri = freeBoundary(TR);
% figure; 
% trisurf(FBtri, points(:,1), points(:,2), points(:,3), ...
%        'FaceColor','cyan','FaceAlpha', 0.8);axis equal;
% TS=triangulation(FBtri, points);
% 
% % Test edge Attachments
% k=1;
% Et = edges(TR);
% Es = edges(TS);
% d=DistanceMatrix(Es, Et); 
% [dm, id]=min(d, [], 2); 
% t2e=id(dm==0); 
% id=true(size(Et,1), 1);
% id(t2e)=false;
% E=Et(id,:);
% ti = edgeAttachments(TR, E);
% figure; hold on;
% tetramesh(tnrb.delaunay, tnrb.points);
% plot3(points(E(k,:),1), points(E(k,:),2), points(E(k,:),3), 'r', 'LineWidth', 2);
% axis equal; view(3);
% 
% figure; hold on;
% tetramesh(tnrb.delaunay(ti{k},:), tnrb.points);
% plot3(points(E(k,:),1), points(E(k,:),2), points(E(k,:),3), 'r', 'LineWidth', 2);
% axis equal; view(3);







