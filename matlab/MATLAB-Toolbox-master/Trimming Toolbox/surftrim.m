function [pts1,pnts1,pts2,pnts2]=surftrim(srf1,srf2,depth,tola)

% Surface trimmed based on equal-depth segment, whthout constructing the
% four-branch tree.

% Input:
%   srf1: NURBS/Bezier surface 1.
%   srf1: NURBS/Bezier surface 1.
%   depth: Current depth of the segment.
%   tola: Pre-depth of the equal-depth segment.
% Output:
%   pts1: Parametric coords of the intersection points. The points are
%       sorted, if there are several lines, each line is stored as a cell
%       array.
%   pnts1: Physical coords of the intersection points.
%   pts2, pnts2: Intersection points of surface 2.
%   d: Distance of the corresponding pair of intersection points, which is
%       used to check the precision.

% Notice: 
% 1.If the 2 box are NOT intersected, then all the outpu variances
% equal [].
% 2.The parametric coords pts1,pts2 are just the evaluated value, which
% are used as the initial value to calculate the precise value.
% 3.The input 2 surfaces must be completed surfaces, meaning they are both
% structure array, without segment into a cell array.

global INTER1 INTER2
% INTER{i}.points=[numpts,3]. INTER{i}.knots=[numpts,2].
% Notice: aviod to repeat to construct box of surface2 and refine surface1.
aabb2=[];
obb2=[];
nbzr1=[];
nbzr2=[];
if depth==1
    INTER1=[];
    INTER2=[];
    
%     flag1=0;
%     flag2=0;
    
    
    % Check if the 2 original input surfaces are intersected using OBB.
    boundbox=nrbboundary(srf1,'whole');
    [vert,~,sa]=OBB_3D('boundary',boundbox,'min');
    box1=boxconstruct_3D('OBB',vert,sa);
    boundbox=nrbboundary(srf2,'whole');
    [vert,~,sa]=OBB_3D('boundary',boundbox,'min');
    box2=boxconstruct_3D('OBB',vert,sa);
    logic=colldetect_3D(box1,box2);
    if logic==0
        pts1=[];pts2=[];pnts1=[];pnts2=[];
        return;
    end
    
    if strcmp(srf1.form,'B-NURBS')
        srf1_=nrb2bzr(srf1);% Cell array.
    elseif strcmp(srf1.form,'Bezier')
        % If the first input surface is Bezier surface.
        srf1_{1}=srf1;
    end
    if strcmp(srf2.form,'B-NURBS')
        srf2_=nrb2bzr(srf2); 
    elseif strcmp(srf2.form,'Bezier')
        % If the first input surface is Bezier surface.   
        srf2_{1}=srf2;
    end
else
    srf1_=srf1;
    srf2_=srf2;
end

num1=numel(srf1_);
num2=numel(srf2_);
for i=1:num1
    % Collision check using AABB
    [ptmin,ptmax]=AABB_3D('Vertex',srf1_{i});
    aabb1=boxconstruct_3D('AABB',ptmin,ptmax);
    for j=1:num2
        % Avoid to construct repeated box.
        if numel(aabb2)<j || isempty(aabb2{j}) 
            [ptmin,ptmax]=AABB_3D('Vertex',srf2_{j});
            aabb2{j}=boxconstruct_3D('AABB',ptmin,ptmax); 
        end
        logic=colldetect_3D(aabb1,aabb2{j});
        if logic==1
            % Collision check using OBB
            boundbox=nrbboundary(srf1_{i});
            [vert,~,sa]=OBB_3D('boundary',boundbox,'min');
            obb1=boxconstruct_3D('OBB',vert,sa);
            % Avoid to construct repeated box.
            if numel(obb2)<j || isempty(obb2{j})
                boundbox=nrbboundary(srf2_{j});
                [vert,~,sa]=OBB_3D('boundary',boundbox,'min');
                obb2{j}=boxconstruct_3D('OBB',vert,sa);  
            end
            logic=colldetect_3D(obb1,obb2{j});
            if logic==1           
                if depth<tola% Need to segment further
                    h1=@(x)((x(1)+x(end))/2);
                    % Avoid to refine the surface repeatedly.
                    if numel(nbzr1)<i || isempty(nbzr1{i})
                        uv=cellfun(h1,srf1_{i}.knots,'UniformOutput',false);                                              
                        nbzr1{i}=bzrseg(srf1_{i},uv);
                    end     
                    if numel(nbzr2)<j || isempty(nbzr2{j})
                        uv=cellfun(h1,srf2_{j}.knots,'UniformOutput',false);
                        nbzr2{j}=bzrseg(srf2_{j},uv);
                    end   
                    % 容易出错：递归调用
                    % Iterative calling to set INTER1, INTER2
                    surftrim(nbzr1{i},nbzr2{j},depth+1,tola);
                else
                    % 2 Patches of the leaf node
                    vert1=srf1_{i}.coefs;
                    vert1=[vert1(:,1,1),vert1(:,end,1),vert1(:,1,end)...
                        vert1(:,end,end)];
                    vert1=vert1(1:3,:)./repmat(vert1(4,:),3,1);
                    vert1=vert1';% Size(vert1)=[4,3]
                    knots1=srf1_{i}.knots;% Cell array
                    vert2=srf2_{j}.coefs;
                    vert2=[vert2(:,1,1),vert2(:,end,1),vert2(:,1,end)...
                        vert2(:,end,end)];
                    vert2=vert2(1:3,:)./repmat(vert2(4,:),3,1);
                    vert2=vert2';
                    knots2=srf2_{j}.knots;  
                    % Whether the 4 triangles are intersected
                    for ii=1:2
                        for jj=1:2
                            % 两个求交曲面地位是对等的
                            [interpt1, interpt2, uv]=intertri2tri(...
                                vert1(ii:ii+2,:),vert2(jj:jj+2,:));
                             % The 2 triangles are coincide      
                            if size(interpt1,1)==4 || size(interpt1,1)==3
%                                 flag1=flag1+1;
                                interpt1=[interpt1(1,:);
                                    sum(interpt1(2:end,:))/(size(interpt1,1)-1)];
                                interpt2=[interpt2(1,:);
                                    sum(interpt2(2:end,:))/(size(interpt2,1)-1)];        
                                % NOTICE: uv is NOT precise.
                                uv=[0.5,0.5;
                                    0.5,0.5;
                                    0.5,0.5;
                                    0.5,0.5;];
                            end
                            if size(interpt1,1)==2                                
                                % Transform global par-coords to local coords.
                                tem1=[];tem2=[];
                                h2=@(u,a,b)((1-u)*a+u*b);
                                for kk=1:2
                                    tem1=[tem1,h2(uv([1,2],kk),...
                                        knots1{kk}(1),knots1{kk}(end))];
                                    tem2=[tem2,h2(uv([3,4],kk),...
                                        knots2{kk}(1),knots2{kk}(end))];
                                end
                                uv=[tem1;tem2];
                                [INTER1,INTER2]=lineconnect(INTER1,INTER2,...
                                    interpt1,interpt2,uv([1,2],:),uv([3,4],:));            
                            elseif size(interpt1,1)==1
%                                 flag2=flag2+1;
                                
                                
                            end       
                        end     
                    end     
                end  
            end    
        end       
    end
    
end

if nargout>1
    h3=@(x)(x.points);
    h4=@(x)(x.knots);                
    tpts1=cellfun(h4,INTER1,'UniformOutput',false);               
    tpnts1=cellfun(h3,INTER1,'UniformOutput',false); 
    tpts2=cellfun(h4,INTER2,'UniformOutput',false); 
    tpnts2=cellfun(h3,INTER2,'UniformOutput',false); 
    while true
        % Connect the segments which can be combined into one line.
        [pnts1,pnts2,pts1,pts2]=conncetline(tpnts1,tpnts2,tpts1,tpts2);
        if numel(pnts1)==numel(tpnts1)
            break;
        else
            tpnts1=pnts1;
            tpnts2=pnts2;
            tpts1=pts1;
            tpts2=pts2;
        end 
    end
end

end

function [ninter1,ninter2]=lineconnect(inter1,inter2,linep1,linep2,linek1,linek2)

% Check if the new line segment belongs to the lines stored in inter( inter
% is a cell array). If it belongs, the line is connected with inter{i} in
% the output ninter, if not,  the 2 points of the line are stored as a new 
% line in ninter.
num=numel(inter1);
T=1e-8;
ninter1=inter1;
ninter2=inter2;
for i=1:num
    pt1=inter1{i}.points(1,:);
    pt2=inter1{i}.points(end,:);
    k1=inter1{i}.knots(1,:);
    k2=inter1{i}.knots(end,:);     
    pt1_=inter2{i}.points(1,:);
    pt2_=inter2{i}.points(end,:);
    k1_=inter2{i}.knots(1,:);
    k2_=inter2{i}.knots(end,:); 
    % Topology of inter2 depends on inter1
    for j=1:2
        if norm(pt1-linep1(j,:))<T
            ninter1{i}.points(1,:)=(pt1+linep1(j,:))/2;
            ninter1{i}.knots(1,:)=(k1+linek1(j,:))/2;
            ninter1{i}.points=[linep1(abs(j-3),:);ninter1{i}.points];
            ninter1{i}.knots=[linek1(abs(j-3),:);ninter1{i}.knots];
            ninter2{i}.points(1,:)=(pt1_+linep2(j,:))/2;
            ninter2{i}.knots(1,:)=(k1_+linek2(j,:))/2;
            ninter2{i}.points=[linep2(abs(j-3),:);ninter2{i}.points];
            ninter2{i}.knots=[linek2(abs(j-3),:);ninter2{i}.knots];
            return;
        elseif norm(pt2-linep1(j,:))<T
            ninter1{i}.points(end,:)=(pt2+linep1(j,:))/2;
            ninter1{i}.knots(end,:)=(k2+linek1(j,:))/2;
            ninter1{i}.points=[ninter1{i}.points;linep1(abs(j-3),:)];
            ninter1{i}.knots=[ninter1{i}.knots;linek1(abs(j-3),:)];
            ninter2{i}.points(end,:)=(pt2_+linep2(j,:))/2;
            ninter2{i}.knots(end,:)=(k2_+linek2(j,:))/2;
            ninter2{i}.points=[ninter2{i}.points;linep2(abs(j-3),:)];
            ninter2{i}.knots=[ninter2{i}.knots;linek2(abs(j-3),:)];            
            return;         
        end
    end
end
% If the line doesn't belong to the existing segments, then create a new
% segment.
ninter1{end+1}.points=linep1;
ninter1{end}.knots=linek1;
ninter2{end+1}.points=linep2;
ninter2{end}.knots=linek2;
end

function [npnts1,npnts2,npts1,npts2]=conncetline(pnts1,pnts2,pts1,pts2)

% Conncet intersection segments into their corresponding segment as
% possible. The input must be cell array and all of them are in accordance.
T=1e-8;
npnts1{1}=pnts1{1};
npnts2{1}=pnts2{1};
npts1{1}=pts1{1};
npts2{1}=pts2{1};
num1=numel(pnts1);

for i=2:num1% Segment needed to be checked
    p1=pnts1{i}(1,:);% S
    p2=pnts1{i}(end,:);% E
    num2=numel(npnts1);
    for j=1:num2% Segment used as criteria
        np1=npnts1{j}(1,:);% S_
        np2=npnts1{j}(end,:);% E_
        if norm(p1-np1)<T% E->S/S_->E, E:end point. S:starting point.
            npnts1{j}=[pnts1{i}(end:2,:);npnts1{j}];
            npnts2{j}=[pnts2{i}(end:2,:);npnts2{j}];
            npts1{j}=[pts1{i}(end:2,:);npts1{j}];
            npts2{j}=[pts2{i}(end:2,:);npts2{j}];
            break;
        elseif norm(p1-np2)<T% S_->E_/S->E
            npnts1{j}=[npnts1{j};pnts1{i}(2:end,:)];
            npnts2{j}=[npnts2{j};pnts2{i}(2:end,:)];
            npts1{j}=[npts1{j};pts1{i}(2:end,:)];
            npts2{j}=[npts2{j};pts2{i}(2:end,:)];
            break;
        elseif norm(p2-np1)<T% S->E/S_->E_
            npnts1{j}=[pnts1{i};npnts1{j}(2:end,:)];
            npnts2{j}=[pnts2{i};npnts2{j}(2:end,:)];
            npts1{j}=[pts1{i};npts1{j}(2:end,:)];
            npts2{j}=[pts2{i};npts2{j}(2:end,:)];                    
            break;  
        elseif norm(p2-np2)<T% S_->E_/E->S
            npnts1{j}=[npnts1{j};pnts1{i}(end-1:1,:)];
            npnts2{j}=[npnts2{j};pnts2{i}(end-1:1,:)];
            npts1{j}=[npts1{j};pts1{i}(end-1:1,:)];
            npts2{j}=[npts2{j};pts2{i}(end-1:1,:)];            
            break;
        end
        if j==num2% Store the checking segment as a new segment in np.
            npnts1{end+1}=pnts1{i};
            npnts2{end+1}=pnts2{i};
            npts1{end+1}=pts1{i};
            npts2{end+1}=pts2{i};     
        end
    end    
end

end

%% demo
% clear;clear global;
% srf1 = nrbcylind(10,8,[],3*pi/2,pi);
% srf2=nrbtestsrf;
% depth=1;
% tola=4;
% [pts1,pnts1,pts2,pnts2]=surftrim(srf1,srf2,depth,tola);
% 
% %%
% figure;hold on;axis equal;view(3);
% nrbplot(srf1,[30,30],'light','on');
% nrbplot(srf2,[30,30],'light','on');
% 
% %%
% for i=1:length(pnts1)
%     ra=rand(1,3);
%     plot3(pnts1{i}(:,1),pnts1{i}(:,2),pnts1{i}(:,3),'*-','color',ra);
% end


% figure;hold on;
% for i=1:numel(bzrsrf)
%    
%     h(i)=nrbplot(bzrsrf{i},[50,50],'light','on');
%   
% end
% for j=1:numel(h)
%     c=rand(1,3);
%     set(h(j),'FaceColor',c);
% end
% axis equal;

