function [ cpolyg,idtemcoefs,nknt,nrb ] = nrbcrvbox_2D( nrb,varargin)

% Construct 2D box based on the control polygon of the NURBS curve, or the
% special points, which depends on the 'option'.

% Input:
%     nrb: The NURBS structure of a curve.
%     nbox: Number of the boxes which are to be constructed. The default
%       number equals the numbers of inertial break points +1.
%     option: 'Control' for constructing boxes based on the control polygon of NURBS, 
%         'Interval' for constructing boxes whose edges are parallel to axes.
% Output:
%   cpolyg: A cell array, which stores control points of each polygon. 
%       The control points have been sorted and the first point is the same
%       as the last one, because cpolyg=coefs(idtemcoefs,:).
%       If the control polygon is a AABB rectangle obtained by setting option as 'Interval',then just stores [Xmin,Ymin] and [Xmax,Ymax] like AABB_2D. 
%       And the number of elements equals the number of polygons. In each cell, number of rows equals number of control points.
%   idtemcoefs: Cell array. Indices of each box's vertices. cpolyg=coefs(idtemcoefs,:), see boxknt.m.
%   nknt: A cell array which stores knots of each part of the curve. If
%       the curve is inserted new knots, then nknt stores new knot segments of new curve's knot vector.
%   nrb: New nurbs structure after the knot insertion.


[knt,order,coefs,brk,ic,rep,stagid]=predeal(nrb); % order= nrb.order-1

if nargin==2
    if ischar(varargin{1})        
        nbox=length(brk)-1;
        option=varargin{1};
    else
        nbox=varargin{1};
        option='Control';
    end
elseif nargin==3
    nbox=varargin{1};
    option=varargin{2};
else
    nbox=length(brk)-1;
    option='Control';
end

insknt=[];
insrep=order-rep;
insid=find(insrep>0);% The break points that are NOT stagnation points.
numinsid=length(insid);
numstagid=length(stagid);  % Number of the break points that are stagnation points.
remstag=nbox-numstagid-1;

if strcmp(option,'Control')
    if  remstag>0
        % Need to insert new knots.       
        if  remstag<=numinsid 
            % Only use the break points
            for i=1:remstag
                insknt=[insknt,repmat(brk(insid(i)),[1,insrep(insid(i))])]; % Knots which need to be insterted.
            end   
        else
            % Insert new knots besides break points.
            restins=remstag-numinsid;
            tembrk=brk;
            for i=1:numinsid
                insknt=[insknt,repmat(brk(insid(i)),[1,insrep(insid(i))])]; % Knots which need to be insterted.
            end   
            while restins>0 % Insert knots by adaptive bisection of each interval.
                for j=1:length(tembrk)-1
                    newknt=(brk(j)+brk(j+1))/2;
                    tembrk=sort([tembrk,newknt]);
                    insknt=[insknt,repmat(newknt,[1,order])]; % Knots needed to be inserted.
                    insknt=sort(insknt);
                    restins=restins-1;
                    if restins<=0
                        break;
                    end
                end
                brk=tembrk;
            end
        end
        nrb = nrbkntins(nrb,insknt);
        [nknt,idtemcoefs,cpolyg]=boxknt(nrb);

    else % The number of boxes is no more than the Bezier segments (which means nbox<numstagid+1).
        % There's no need to insert new knots, but just use the stagnations.
        [nknt,idtemcoefs,cpolyg]=boxknt(nrb,nbox);    
    end
        
elseif strcmp(option,'Interval')
    % Use the special points, whose slopes are 0 or infinite or are stagnation points, to construct
    % rectanglar boxes. Insert knots and control ponits when refining the boxes if necessary.
    stag=brk(stagid);
    pts=nrbcrvslope(nrb);
    pts=unique([pts,0,1,stag]); % All special points' par cords.
    numpts=length(pts);
   
    if nbox>numpts-1
    % Need to insert new knots.
    % Divide each interval based on the existed boxes adaptively.
        restins=nbox-numpts+1;
        tempts=pts; % Break points.
        while restins>0
            for j=1:length(tempts)-1
                newknt=(pts(j)+pts(j+1))/2;
                % If newknt is the break point which has already
                % existed, then the times of insertion is order-rep.
                teminsid=find(abs(brk-newknt)<1e-6);
                repeat=rep(teminsid);
                if isempty(repeat)
                    tempts=sort([tempts,newknt]);
                    insknt=[insknt,repmat(newknt,[1,order])];
                    insknt=sort(insknt);
                else
                    insknt=[insknt,repmat(newknt,[1,order-repeat])];
                    insknt=sort(insknt);
                end
                restins=restins-1;
                if restins<=0
                    break;
                end
            end
            pts=tempts;
        end
    end    
    nrb = nrbkntins(nrb,insknt); % If insknt=[], nrb doesn't change. 
    knts=nrb.knots;
    pnts=nrbeval(nrb,pts); 
    pnts(3,:)=[];pnts=pnts'; % [numpnts,2]
    startid=1;
    for i=1:nbox-1
        [p1,p2,~,idtemcoefs{i}]=AABB_2D(pnts(i:i+1,:));
        cpolyg{i}=[p1(:),p2(:)]'; % [Xmin Xmax;Ymin Ymax]
        temid=find(knts==pts(i+1));
        nknt{i}=knts(startid:temid+order-1);
        startid=temid;
    end
        [p1,p2,~,idtemcoefs{i}]=AABB_2D(pnts(nbox:end,:));
        cpolyg{nbox}=[p1(:),p2(:)]'; % [Xmin Xmax;Ymin Ymax]
        nknt{i}=knts(startid:end);    
else
    warning('The method of constructing box is wrong');
end

end

function [nknt,idtemcoefs,cpolyg]=boxknt(varargin) % Use the function convhull.m
 % Get the new knots and control polygons (repreented by their coordinates
 % of polygons' vertices) from the inputting nrb and number of box. If
 % there's the second input nbox, nbox< numstagid, otherwise just input one
 % argument without nbox.
        nrb=varargin{1};
        [knt,order,coefs,brk,ic,rep,stagid]=predeal(nrb);     
        numstagid=length(stagid);  
        nbox=numstagid+1;
        begid=1;
        idp1=1;
        if nargin==2
            nbox=varargin{2};
        end      
        for i=1:nbox-1
            endid=ic(stagid(i))+rep(stagid(i))-1;
            nknt{i}=knt(begid:endid);
            begid=ic(stagid(i));          
            p=nrbeval(nrb,brk(stagid(i)));p(3)=[];p=p';
            idp2=find(abs(coefs(:,1)-p(1))<10e-6 & abs(coefs(:,2)-p(2))<10e-6);
            temcoefs=coefs(idp1:idp2,:);
            idtemcoefs{i}=convhull(temcoefs(:,1),temcoefs(:,2));
            % The output is a dim*1 vector, each element represents the index of a point. The first element is the same as the last one.
            cpolyg{i}=temcoefs(idtemcoefs{i},:);
            idp1=idp2;
        end
        nknt{end+1}=knt(begid:end);      
        temcoefs=coefs(idp1:end,:);
        idtemcoefs{end+1}=convhull(temcoefs(:,1),temcoefs(:,2));
        cpolyg{end+1}=temcoefs(idtemcoefs{end},:);%cpolyg{end+1}=temcoefs(idtemcoefs{end+1},:);
    end

function [knt,order,coefs,brk,ic,rep,stagid]=predeal(nrb)
% Inputing nrb, get the knots, order, physical coordinates of control points (represented
% by [numpoints,2] vector), break points and their indices and their corresponding
% repeatabilities, indices of stagnation points.
    knt=nrb.knots;
    order=nrb.order-1;
    coefs=nrb.coefs;
    % Transform the homogeneous cords to card cords.
    coefs=coefs(1:3,:)./repmat(coefs(4,:),[3,1]);
    coefs=coefs([1,2],:);coefs=coefs'; % Number of rows equals number of points.
    % Get all the break points and their corresponding repeaties.
    [brk,ic]=unique(knt);
    rep=[ic(2:end);length(knt)+1];
    rep=rep-ic;
    stagid=find(rep==order); % Stagnation points' indices corresponding to rep/brk/ic, except two ending points.
end

%% demo
% crv1=nrbtestcrv;
% crv2=nrbcirc;
% coefs =[ 6.0  0.0  6.0  1;
%          -5.5  0.5  5.5  1;
%          -5.0  1.0 -5.0  1;
%           4.5  1.5 -4.5  1;
%           4.0  2.0  4.0  1;
%          -3.5  2.5  3.5  1;
%          -3.0  3.0 -3.0  1;
%           2.5  3.5 -2.5  1;
%           2.0  4.0  2.0  1;
%          -1.5  4.5  1.5  1;
%          -1.0  5.0 -1.0  1;
%           0.5  5.5 -0.5  1;
%          0.0  6.0  0.0  1]';
%  knots = [0 0 0 0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1 1 1 1];
%  coefs(3,:)=0;
%  crv3 = nrbmak(coefs,knots);
%  
%  [c1,n1]=nrbcrvbox_2D(crv1);
%  figure;
%  nrbplot(crv1,100);
%  hold on;
%  num=length(c1);
%  for i=1:num
%      plot(c1{i}(:,1),c1{i}(:,2));
%  end
%  
%  [c2,n2]=nrbcrvbox_2D(crv2);
%  figure;
%  nrbplot(crv2,100);
%  hold on;
%  num=length(c2);
%  for i=1:num
%      plot(c2{i}(:,1),c2{i}(:,2));
%  end
%  
%  [c3,n3]=nrbcrvbox_2D(crv3);
%  figure;
%  nrbplot(crv3,100);
%  hold on;
%  num=length(c3);
%  for i=1:num
%      plot(c3{i}(:,1),c3{i}(:,2));
%  end
%  
%   [c1,n1]=nrbcrvbox_2D(crv1,'Interval');
%  figure;
%  nrbplot(crv1,100);
%  hold on;
%  num=length(c1);
%  for i=1:num
%      pt1=c1{i}(1,:);pt2=c1{i}(2,:);
%      pt=[pt1;pt2(1) pt1(2);pt2;pt1(1) pt2(2);pt1];
%      plot(pt(:,1),pt(:,2));
%  end 
%  
%  [c2,n2]=nrbcrvbox_2D(crv2,'Interval');
%  figure;
%  nrbplot(crv2,100);
%  hold on;
%  num=length(c2);
%  for i=1:num
%      pt1=c2{i}(1,:);pt2=c2{i}(2,:);
%      pt=[pt1;pt2(1) pt1(2);pt2;pt1(1) pt2(2);pt1];
%      plot(pt(:,1),pt(:,2));
%  end 
%
%  [c3,n3]=nrbcrvbox_2D(crv3,'Interval');
%  figure;
%  nrbplot(crv3,100);
%  hold on;
%  num=length(c3);
%  for i=1:num
%      pt1=c3{i}(1,:);pt2=c3{i}(2,:);
%      pt=[pt1;pt2(1) pt1(2);pt2;pt1(1) pt2(2);pt1];
%      plot(pt(:,1),pt(:,2));
%  end 
