% Transform the directions of plane surfaces to have a common normal vector.
% 
%   Input:  
%       srfs - a struct array of nurbs surfaces with arbitrary normal direction
%       dr -  a refference direction (default is [0 0 1])
%        
%   Output :  
%      srfs - a struct array of nurbs surfaces with the same normal direction

function srfs=nrlsrfsdirection(srfs, dr)

if nargin==1
    dr= [0 0 1];
end

% Apply the refference direction
[~, ~, ~, dr1] = extractsrf(srfs(1));
dn=dot(dr1(:,1,1)', dr);
if dn<0
    srfs(1)=nrltransp(srfs(1));
end

% Get a struct of the points and directions of the surfaces
m=numel(srfs);
srfeds=cell(4,m);
for k=1:m
    [srfeds{1,k}, srfeds{2,k}, srfeds{3,k}, srfeds{4,k}] = extractsrf(srfs(k));
end

% Get all points and the corresponding surface directions
pnts=zeros(4*m,3); dpns=pnts;
for i=1:m
    pnts(4*i-3,:)=srfeds{1,i}(:,1,1);
    pnts(4*i-2,:)=srfeds{1,i}(:,2,1);
    pnts(4*i-1,:)=srfeds{1,i}(:,1,2);
    pnts(4*i-0,:)=srfeds{1,i}(:,2,2);
    dpns(4*i-3,:)=srfeds{4,i}(:,1,1);
    dpns(4*i-2,:)=srfeds{4,i}(:,2,1);
    dpns(4*i-1,:)=srfeds{4,i}(:,1,2);
    dpns(4*i-0,:)=srfeds{4,i}(:,2,2);
end

% Correct the direction of surfaces
flag=zeros(m,1); flag(1)=1;
flag=flag==1;
dm = DistanceMatrix(pnts, pnts);
tol=(norm(pnts(1,:)-pnts(2,:)))*1e-3;
dm=(dm<tol); 
for i=1:m
    for j=1:m
        for s=1:4
            for t=1:4
                p=4*i-s+1; q=4*j-t+1; 
                if p~=q && dm(p,q)
                    dn=dot(dpns(p,:), dpns(q,:));
                    if dn<0
                        if flag(i) && ~flag(j)
                            srfs(j)=nrltransp(srfs(j));
                            [srfeds, dpns]=updatesrfeds(srfs, srfeds, dpns, j);
                            flag(j)=true;
                        elseif ~flag(i) && flag(j)
                            srfs(i)=nrltransp(srfs(i));
                            [srfeds, dpns]=updatesrfeds(srfs, srfeds, dpns, i);
                            flag(i)=true;
                        end
                    elseif dn>0
                        if flag(i) && ~flag(j)
                            flag(j)=true; 
                        elseif ~flag(i) && flag(j)
                            flag(i)=true; 
                        end
                    end
                end
            end
        end
    end
end

% Extract the four coner points and their derivatives of a surface
function varargout =extractsrf(varargin)

srf=varargin{1};
ut = [0, 1]; vt = [0, 1]; 
[pnts, jac] = nrldeval(srf,{ut,vt}); 
varargout{1} = pnts; 
varargout{2} = jac{1}; 
varargout{3} = jac{2}; 
varargout{4} = cross(jac{1}, jac{2}, 1);

% Update the above information of a surface
function [srfeds, dpns]=updatesrfeds(srfs, srfeds, dpns, i)

[srfeds{1,i}, srfeds{2,i}, srfeds{3,i}, srfeds{4,i}] = extractsrf(srfs(i));
dpns(4*i-3,:)=srfeds{4,i}(:,1,1);
dpns(4*i-2,:)=srfeds{4,i}(:,2,1);
dpns(4*i-1,:)=srfeds{4,i}(:,1,2);
dpns(4*i-0,:)=srfeds{4,i}(:,2,2);

