function bzrsegment(mode,varargin)

% Segment the NURBS/Bezier model by 2 modes: 'Adaptive' for adaptive segment and
% 'Equal' for equal-depth segment. Nocitce that this function is only used
% to segment the surface and get the tree structure, without trimming or
% mesh generation.

% Input in old edition: model,tol,mode,depth,varargin
% Input in new edition: mode,varargin. For 'adaptive': model=varargin{1},
% tol=varargin{2}, depth =varargin{3}, parentid_=varargin{4}, currentid_=
% varargin{5}; For 'equal': model1=varargin{1}, model2=varargin{2},
% depth=varargin{3}, tola=varargin{4}, parentid_=
% varargin{5}, currentid_=varargin{6}

% Input:
%   model: NURBS/Bezier model, including curve or surface, representing one node 
%       in the tree (adaptive) or the 4 NODES in one depth (for equal-depth),
%       which is a cell array.And each cell stands a Bezier model.
%   mode:'Adaptive' or 'Equal'.
%   depth: Depth of the current iteration.
%   tol: Ttolerance for adaptive segment of checking the flatness,
%       see flatcheck.m.
%   tola: Limit of segment depth for 'equal' mode.
%   parentid_: Parent node index of the current node.
%   currentid_: Index if the current node(the input node needed to be 
%       processed), used to record the 1st node of each depth.
% Output:
%   modeltree1: Multiway tree of the whole model,which is a cell array for
%       'adaptive' mode.
%   modeltree1, modeltree2: Multiway tree of the 2 original input surfaces,
%       which is 2 cell arrays for  'equal' mode.

% modeltree1.model/node: Use the 1*4 vector to represent the tree structure.
% modeltree1 is a cell array, the number of modeltree1 equals the number of
% all the nodes of the tree.

global NUMNODE modeltree 
% clear global; clear global var;
if strcmp(mode,'Adaptive')   
    model=varargin{1};
    tol=varargin{2};
    depth=varargin{3};   
    if nargin==5
        parentid_=varargin{4};
    elseif nargin==6
        parentid_=varargin{4}; 
        currentid_=varargin{5}; 
    end
elseif strcmp(mode,'Equal')
    model1=varargin{1};
    model2=varargin{2};
    depth=varargin{3};
    tola=varargin{4};
    if nargin==6
        parentid_=varargin{5};
    elseif nargin==7
        parentid_=varargin{5}; 
        currentid_=varargin{6}; 
    end
else
    error('The mode must be ''Adaptive'' or ''Equal''');
end
           
if depth==1
    % Prepare for another tree.
    NUMNODE=[];
    modeltree=[];
    if strcmp(mode,'Adaptive')
        if iscell(model) 
    % The original input is a multiple cell array,eg. several patches after
    % segment a NURBS surface
            depth=2;
            [num1,num2]=size(model);
            if isempty(NUMNODE)
                NUMNODE=0;
            else
                modeltree{NUMNODE}.node(4)=num1*num2;
            end        
            currentid_=NUMNODE+1;% Index of 1st node after further segment.
            parentid_=1;
            for i=1:num1
                for j=1:num2
                    modeltree{NUMNODE+1}.node=[depth,1,i*10+j];
                    modeltree{NUMNODE+1}.model=model{i,j};
                    NUMNODE=NUMNODE+1;
                end
            end
        else 
        % The original input is a structure, meaning it's only one patch
        % (Bezier or NURBS).       
            NUMNODE=1;       
            modeltree{NUMNODE}.model=model;
            modeltree{NUMNODE}.node=[depth,0,11];
            if strcmp(model.form,'B-NURBS')
                % If the original input model is NURBS model, 
                % then transform it into Bezier model.      
                nmodel=nrb2bzr(model);         
                bzrsegment('Adaptive',nmodel,tol,depth,1,2);
                return;
            end
        end
    else% 'Equal'
        % If mode='equal', default condition is that the 2 original input 
        % surfaces are structure array of NURBS/Bezier model. If there's 
        % at least one NURBS model, then both surfaces' depth=2.
        NUMNODE=1;       
        modeltree{NUMNODE}.model=model1;
        modeltree{NUMNODE}.node=[depth,0,1];
        if strcmp(model1.form,'B-NURBS')   
            nmodel1=nrb2bzr(model1);    
        end
        if strcmp(model2.form,'B-NURBS')   
            nmodel2=nrb2bzr(model2);    
        end        
        num1=numel(nmodel1);
        modeltree{NUMNODE}.node(4)=num1;       
        for i=1:num1
            modeltree{NUMNODE+1}.node=[depth+1,1,i];
            modeltree{NUMNODE+1}.model=nmodel1{i};
            NUMNODE=NUMNODE+1;
        end
        
        % 递归
        bzrsegment('Equal',nmodel1,nmodel2,depth+1,tola,1,2);
        return;
    end               
end
        
if strcmp(mode,'Adaptive')       
    [num1,num2]=size(model);% num1: Number of ros. num2: Number of columns.
    if num1*num2==1 % If there's one structure, transform to cell array.
        model_{1,1}=model;
    else
        model_=model;
    end    
    % The parent index is also represented by vector-form. Each parent 
    % node is the current node of the last depth.
    currentid=currentid_;    
    for ii=1:num1
        for jj=1:num2
            [uv,lg]=flatcheck(model_{ii,jj},tol);      
            if lg==1
                % nmodel: Child node of the current parent node. 
                nmodel=bzrseg(model_{ii,jj},uv);               
                [num3,num4]=size(nmodel);
                % 当前节点编号采用矢量表示形式，只适应于曲面 
                % Every node stores a Bezier patch, meaning
                % modeltree1{i}.model is a structure.      
                temcurrentid=NUMNODE+1;                 
                for mm=1:num3
                    for nn=1:num4
                        modeltree{NUMNODE+1}.node=[depth+1,currentid,...
                            mm*10+nn];
                        modeltree{NUMNODE+1}.model=nmodel{mm,nn};
                        NUMNODE=NUMNODE+1;  
                    end
                end
                % Add the number of child nodes in each current node: the
                % 4th element in the tree.node.               
                id1=findmodeltree(modeltree,depth,'depth');
                id3=findmodeltree(modeltree,parentid_,'parent');
                id2=findmodeltree(modeltree,ii*10+jj,'current');
                id=intersect(id1,id2);id=intersect(id,id3); 
                for numid=1:length(id)
                    modeltree{id(numid)}.node(4)=num3*num4;
                end                   
                % 递归调用:针对每个当前结点（父节点）所对应的曲面片
                bzrsegment('Adaptive',nmodel,tol,depth+1,currentid_,...
                    temcurrentid);               
                currentid=currentid+1;
                currentid_=currentid_+1;  
            else
            % The current node is leaf node.Add the number of child 
            % nodes in each current node: the 4th element in the tree.node.
                id1=findmodeltree(modeltree,depth,'depth');
                id3=findmodeltree(modeltree,parentid_,'parent');
                id2=findmodeltree(modeltree,ii*10+jj,'current');
            % parentid_: Index of parent node;
            % currentid: Index of the current node.
                id=intersect(id1,id2);id=intersect(id,id3);  
                for numid=1:length(id)
                    modeltree{id(numid)}.node(4)=0;
                end                   
                currentid=currentid+1;  
                currentid_=currentid_+1;                  
            end
        end              
    end       
   
elseif strcmp(mode,'Equal')
    aabb2=[];
    obb2=[];
    nbzr1=[];
    nbzr2=[];   
    num1=numel(model1);
    if num1==1 
        model_1{1}=model1;
    else
        model_1=model1;
    end         
    num2=numel(model2);
    if num2==1 
        model_2{1}=model2;
    else
        model_2=model2;
    end       
    if depth<tola
        for i=1:num1
            tem1=0;
            [ptmin,ptmax]=AABB_3D('Vertex',model_1{i});
            aabb1=boxconstruct_3D('AABB',ptmin,ptmax);
            for j=1:num2
                if numel(aabb2)<j || isempty(aabb2{j}) 
                    [ptmin,ptmax]=AABB_3D('Vertex',model_2{j});
                    aabb2{j}=boxconstruct_3D('AABB',ptmin,ptmax); 
                end
                logic=colldetect_3D(aabb1,aabb2{j});                
                if logic==1 % AABB-AABB intersected
                    boundbox=nrbboundary(model_1{i});
                    [vert,~,sa]=OBB_3D('boundary',boundbox,'min');
                    obb1=boxconstruct_3D('OBB',vert,sa);
                    if numel(obb2)<j || isempty(obb2{j})
                        boundbox=nrbboundary(model_2{j});
                        [vert,~,sa]=OBB_3D('boundary',boundbox,'min');
                        obb2{j}=boxconstruct_3D('OBB',vert,sa);  
                    end
                    logic=colldetect_3D(obb1,obb2{j});                  
                    if logic==1 % OBB-OBB intersected
                        h1=@(x)((x(1)+x(end))/2);
                        temcurrentid=NUMNODE+1;                        
                        if numel(nbzr2)<j || isempty(nbzr2{j})
                            uv2=cellfun(h1,model_2{j}.knots,...
                                'UniformOutput',false);
                            nbzr2{j}=bzrseg(model_2{j},uv2);                                                                         
                        end                      
                        if numel(nbzr1)<i || isempty(nbzr1{i})
                            uv1=cellfun(h1,model_1{i}.knots,...
                                'UniformOutput',false);                                              
                            nbzr1{i}=bzrseg(model_1{i},uv1);
                            num3=numel(nbzr1{i});
                            for mm=1:num3
                                modeltree{NUMNODE+1}.node=...
                                    [depth+1,currentid_,mm];
                                modeltree{NUMNODE+1}.model=nbzr1{i}{mm};
                                NUMNODE=NUMNODE+1; 
                            end
                            id1=findmodeltree(modeltree,depth,'depth');
                            id3=findmodeltree(modeltree,parentid_,'parent');
                            id2=findmodeltree(modeltree,i,'current');
                            id=intersect(id1,id2);id=intersect(id,id3);                        
                            modeltree{id}.node(4)=4; 
                            
                            % 递归调用
                            bzrsegment('Equal',nbzr1{i},nbzr2{j},depth+1,...
                                tola,currentid_,temcurrentid);   
                        end                       
                    else % OBB-OBB NOT intersected
                        if tem1==num2-1
                            id1=findmodeltree(modeltree,depth,'depth');
                            id2=findmodeltree(modeltree,i,'current');
                            id3=findmodeltree(modeltree,parentid_,'parent');
                            id=intersect(id1,id2);id=intersect(id,id3);  
                            if numel(modeltree{id}.node)==3 
                                modeltree{id}.node(4)=0; 
                            end                   
                        else
                            tem1=tem1+1;
                        end
                    end                   
                else % AABB-AABB NOT intersected
                    if tem1==num2-1
                        id1=findmodeltree(modeltree,depth,'depth');
                        id2=findmodeltree(modeltree,i,'current');
                        id3=findmodeltree(modeltree,parentid_,'parent');
                        id=intersect(id1,id2);id=intersect(id,id3); 
                        if numel(modeltree{id}.node)==3 
                            modeltree{id}.node(4)=0; 
                        end
                    else
                        tem1=tem1+1;                
                    end
                end              
            end     
            currentid_=currentid_+1;   
        end       
    else% depth>tola
        for i=1:num1
            id1=findmodeltree(modeltree,depth,'depth');
            id2=findmodeltree(modeltree,i,'current');
            id3=findmodeltree(modeltree,parentid_,'parent');
            id=intersect(id1,id2);id=intersect(id,id3); 
            try
                if numel(modeltree{id}.node)==3 
                    modeltree{id}.node(4)=0; 
                end               
            catch
                keyboard;
            end
        end   
    end
       
else
    error('The mode is neither ''Adaptive'' nor ''Equal''');
end

end

%% demo
% clear global;
% global NUMNODE modeltree
% % bzrsegment('Adaptive',bzr,1,1);
% % bzrsegment('Adaptive',bzrsrf,1,1);
% bzrsegment('Adaptive',srf,1,1);

%% demo
% clear global;
% global NUMNODE modeltree
% srf1=nrbtestsrf;
% srf2 = nrbcylind(10,8,[],3*pi/2,pi);
% 
% bzrsegment('Equal',srf1,srf2,1,7);
% trimodel=segmenttriangulation(modeltree);
% 
% figure;
% trisurf(trimodel.triangulation,trimodel.points(:,1),trimodel.points(:,2),...
%     trimodel.points(:,3));
% axis equal;
% figure;
% hold on;
% knots=trimodel.knots;
% plot(knots(:,1),knots(:,2),'ro');
% triplot(trimodel.triangulation,knots(:,1),knots(:,2));
% axis equal;







