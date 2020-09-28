function tri=moletriangulation(nrb,varargin)

% Triangulation based on the molecular dynamics method. The resulting mesh
% elements are relatively uniform which means almost all the edges of the
% elements are the same. Notice that this method is not a adaptive method.

% Input:
%   nrb: NURBS structure representing a curve or a surface.
%   h0: Expected length of the mesh edges(mesh seed).
%   m,q: Mass and charge of all the points.The default is 1.
%   codis: Coefficient in the denominator which influences the efficient distance
%       of the neighbor points. codis>=1 and the default is codis=2.
%   cocur: Normalization constant coefficient which is used to calculate the repeling force coefficient. The default is 1.
%   codam: Damping coefficient and the default is 2.
%   covel: Order of the velocity to calculate the damping force. The default is 2.
% Output:
%   tri: Triangular mesh structure after DT triangulation.

m=1;
q=0.016;
codis=2;
cocur=1;
codam=100;
covel=2;
Uformer=0;% Initial potential energy of the system
delta_time=0.01;% Time step to calculate the velocity using acceleration
N=20;% Change the default number of the mesh. N means the number in each iso-parametric line.
% de=true;


% Deal with the degraded points (eg. sphere surface): 1.Coincident points'
% force are 0; Or 2. Preprocessing the degraded points(see Uniq-function).
% 非正则曲面的情况：参数点和老师的几乎完全一样，用老师的点进行nrbeval后得到的points画也不对，可能是nrb2tri函数里points的计算方法不是单纯的nrbeval(nrb,knt)

% NURBS surface. Still needs NURBS curve condition.
[Lu, Lv]=nrbsrfmeasure(nrb);
ru=max(Lu)/min(Lu);
rv=max(Lv)/min(Lv);
if nargin==1
    h0=max([Lu, Lv])/N;
elseif nargin==2
    h0=varargin{1};
elseif nargin==3
    h0=varargin{1};
    m=varargin{2};
elseif nargin==4
    h0=varargin{1};
    m=varargin{2};
    q=varargin{3};
elseif nargin==5
    h0=varargin{1};
    m=varargin{2};
    q=varargin{3};
    codis=varargin{4};
elseif nargin==6
    h0=varargin{1};
    m=varargin{2};
    q=varargin{3};
    codis=varargin{4};    
    cocur=varargin{5};
elseif nargin==7
    h0=varargin{1};
    m=varargin{2};
    q=varargin{3};
    codis=varargin{4};    
    cocur=varargin{5};
    codam=varargin{6}; 
elseif nargin==8
    h0=varargin{1};
    m=varargin{2};
    q=varargin{3};
    codis=varargin{4};    
    cocur=varargin{5};
    codam=varargin{6}; 
    covel=varargin{7}; 
end
% Adapt the u and v direction
if ru<rv
    nrb=nrbtransp(nrb);
    Lv=Lu;% Make sure Lv represent the uniform distance isoparameteric line
end
mm=round(max(Lv)/h0); 
nodes=cell(1, mm); 
s=linspace(0, 1, mm);% Knots are distributed uniformly in par-domain, but not in phy-domain.
[Lu, ~]=nrbsrfmeasure(nrb, s);% Lu represents isoparameteric lines that are non-uniform.
nd=0; 
num=zeros(mm, 1);

% Define the function of calculating the potential energy
Epotential=@(C,codis,q,h)(C*q*q/((codis-1)*h^(codis-1))+10);
% Get the par-cord of initial nodes in par-domain. This can avoid the degraded points.
for i=1:mm
    num(i)=max([ceil(Lu(i)/h0), 1]); % Essential!!!!!!  Number of nodes in each isopar-line in v-direction
    nodes{i}=linspace(0, 1, num(i));% v(y) coordinate of points in each isopar-line in v-direction
    nd=nd+num(i);% Total number of all the nodes
end
x=zeros(nd, 1); % Column vector
y=x; 
t=0;
V=zeros(nd,2);% At time t=0, the velocity is v=0
for i=1:mm
    ni=num(i); 
    % Set the u(x) and v(y) coordinates in par-domain
    x(t+1:t+ni)=s(i);
    y(t+1:t+ni)=nodes{i};
    t=t+ni;
end

delta_space=[abs((s(2)-s(1))/100),abs((nodes{2}(1)-nodes{2}(2))/100)];% Coordinate step to calculate the local curvature.
points=[x,y];

% Repeling force: inside the domain; Attraction force: outside the domain
while true
    % Get the resultant force of each point in par-domain, including points in the internal domain and in the boudnary.
    Force=[];% Resultant force of all the points.
    U=0;% Initial potential energy of the whole system for each time step, which needs to be set to 0 at first.
    for i=1:nd
        force=[0,0];% Resultant force of one pointi.
        Uij=0;% Initial potential energy of one pointi.
        pi=points(i,:);
        % Using the average distance to control the number of the checking points
        Pr=points-repmat(pi,nd,1);
        dist=sqrt(sum(Pr.*Pr,2));
        distol=sum(dist)/(nd-1);% Define the efficient distance of interaction points.
        position=find(dist<=distol & dist>eps);% Find the points in the efficient domain except the pointi itself.
        for j=1:length(position)
            pj=points(position(j),:);       
            lij(1,:)=logical(pi==0 | pi==1);
            lij(2,:)=logical(pj==0 | pj==1);           
            % Needn't to consider the degraded points condition.
            % i-internal, j-boundary or internal
            [force1,C]=force_repel_bound(nrb,pi,pj,delta_space,'repel',codis,q,cocur);
            if sum(lij(1,:))~=0% i-boundary
                force1(lij(1,:))=0;% If the pointi is in the boundary(including the corner vertex), the corresponding component coordinate equals 0.
            end
            % Calculate the force, potential energy of one pointi.
            force=force+force1;
            temUij=Epotential(norm(C),codis,q,norm(pj-pi));
            Uij=Uij+temUij;      
        end
        Force=[Force;force];% Resultant repeling force of all the points.
        U=[U,Uij];% Potential energy of all the points.
    end
    % force1 is the repeling force of one ponit i: dim=1*2;force2 is the resistance force of all the points: dim=nd*2
    force2=force_resist(V,codam,covel);% Resistance force   
    % About the damping coeffficident.
%     if de==true
%         cocur=norm(force2(1,:))/norm(force1(1,:));
%         de=false;
%     end
    Force=Force+force2; % Resultant repeling and resistance force of all the force
    % Determine whether continue to iteration, using the total potential energy of the whole system
    if sum(U)>Uformer
        A=Force/m;% Acceleration of all the points.
        V=V+A*delta_time;
        points=points+V*delta_time;       
        % Compulsorily move the points outside the boundary to the boundary. Needn't to consider the 'attract' force if outside the boundary.
        points(points<0)=0;
        points(points>1)=1;        
        % There may be duplicated points which need to be dealt with.
        [points_,ip]=unique(points,'rows');
        nump=length(points);
        while length(points_)~=nump
            temip=logical(1:nump);
            temip(ip)=0;
            index0=points==0;
            index1=points==1;
            points(temip,:)=points(temip,:)-repmat(delta_space/100,sum(temip),1);
            points(index0)=0;
            points(index1)=1;
            [points_,ip]=unique(points,'rows');
        end

        Uformer=sum(U);
    else% If the Energy at present is lower than that of the former, then stop the iteration
        break;
    end    
end

% Triangulation using the optimized-positioned points.
tri=delaunayTriangulation(points);% points represent the par-cord in par-domain.

end

%% demo
% tic;
% nrb=nrbtestsrf;
% tri=moletriangulation(nrb);
% T=toc;
% figure;
% triplot(tri);
% points=nrbeval(nrb,tri.Points);
% figure;
% trisurf(tri.ConnectivityList,points(1,:),points(2,:),points(3,:));
% title(['Time consuming is: ',num2str(T),' s']);
%%
% x = linspace (-3, 3, 40);
% y = linspace (-3, 3, 40);
% [X, Y] = meshgrid (x, y);
% Z = peaks (X, Y);
% srf1 = bspinterpsurf (X, Y, Z, [2 2], 'equally_spaced');
% tic;
% tri=moletriangulation(srf1);
% T=toc;
% figure;
% triplot(tri);
% points=nrbeval(srf1,tri.Points);
% figure;
% trisurf(tri.ConnectivityList,points(1,:),points(2,:),points(3,:));
% title(['Time consuming is: ',num2str(T),' s']);

%%
% center=[5,5,4];
% circ=nrbcirc(4, center, 0, pi);
% srf1=nrbrevolve(circ, center, [1,0,0], 2*pi);
% tic;
% tri=moletriangulation(srf1);
% T=toc;
% figure;
% triplot(tri);
% points=nrbeval(srf1,tri.Points);
% figure;
% trisurf(tri.ConnectivityList,points(1,:),points(2,:),points(3,:));
% title(['Time consuming is: ',num2str(T),' s']);






