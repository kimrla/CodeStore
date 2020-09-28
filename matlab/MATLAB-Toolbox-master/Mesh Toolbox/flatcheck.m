function [uv,lg]=flatcheck(bzr,tol,mode)

% Check the flatness of a Bezier model(curve or surface) by 2 kinds of
% methods. Notice this function is only suitable for Bezier model.

% Input:
%   bzr: A Bezier structure. NOTICE: size(coefs)=[dim,nu,nv].
%   tol: Tolerance to determine whether the model satisfy the flatness criterion.
%   mode: 'control': Use the control boundary box to calculate the 
%       approximated curvature of the Bezier segment(curve) or patch(surface).
%       'vector': Use the tangent and normal vectors to calculate the
%       curvature-> For an arbitrary parametric model without control boundaries.
% Output:
%   uv: Par-coords of the input Bezier model to insert knots. For curve, uv
%       is a vector, For surface, uv is a cell array {[u],[v]}. uv is
%       par-coords in LOCAL system.
%   lg: Logical scalar. 1: The model doesn't satisfy the flatness and uv is
%       NOT empty ( Still needs to segment the patch further); 0: satisfy, uv is empty.

if nargin==2
    mode='control';
end
% 嵌套函数：子函数可以使用父函数的变量，此时当做全局变量来用，反过来父函数
%不能使用子函数中新定义的变量。如果是当作全局函数，则每处修改都会改变变量的值
points=bzr.coefs;%  NOTICE: size(coefs)=[dim,nu,nv].
ord=bzr.order;
knots=bzr.knots;
% Determine if the model need to be segmented and find the inserted points.
if strcmp(mode,'control')
    if length(ord)==2 % surface
        [uv{1},lg1]=finduvins(points,tol);
        tempt=permute(points,[1,3,2]);
        [uv{2},lg2]=finduvins(tempt,tol);
        lg=lg1 | lg2;
        if lg==1
            [uv,~]=loc_gol(bzr,uv,'local');% Transform global domain to local domain.
        end                        
    elseif length(ord)==1 % curve
        points=points(1:3,:)./repmat(points(4,:),3,1);
        [~,numu]=size(points);
        % Deal with each row or column of the control lines.
        tempoints=points';
        L=[tempoints(1,:);tempoints(end,:)];
        for j=2:numu-1
            % Number of columns equal to number of isopar-lines. 
            % Number of rows equal to number of INTERNAL points in each isopar-line.
            [~,~,d(j-1)]=interpoint2line(L,tempoints(j,:));   
        end
        uvins=chordpar(points);% Each row represent a knot vector of a isopar-line.
        idd=find(d>=tol);
        if isempty(idd)
            lg=0;
            uv=[];
        else
            lg=1;
            uv=uvins(idd);
        end
        if lg==1
            [uv,~]=loc_gol(bzr,uv,'local');% Transform global domain to local domain.
        end   
    end
%%  For an arbitrary model, if it has no control boundaries, then 
% use thetangent and normal vectors to approximate the curvature.
elseif strcmp(mode,'vector')
    
    
%%    
else
    error('The mode need to be ''vector'' or ''control''');
end


    function [uv,lg]=finduvins(points,tol)
    % Find par-points needed to be inserted.
        points=points(1:3,:,:)./repmat(points(4,:,:),3,1);
        [~,numu,numv]=size(points);
        % Deal with each row or column of the control lines.
        for i=1:numv
            tempoints=points(:,:,i);tempoints=tempoints';
            L=[tempoints(1,:);tempoints(end,:)];
            for jj=2:numu-1
                % Number of columns equal to number of isopar-lines.
                % Number of rows equal to number of INTERNAL points in each isopar-line.
                [~,~,dd(jj-1,i)]=interpoint2line(L,tempoints(jj,:));   
            end
            uvins(i,:)=chordpar(points(:,:,i));% Each row represent a knot vector of a isopar-line.
        end
        aved=sum(dd,2)/numv;
        idd=find(aved>=tol);
        if isempty(idd)
            lg=0;
            uv=[];
        else
            lg=1;
            uvins=sum(uvins)/numv;
            uv=uvins(idd);
        end
    end

    function uv=chordpar(A)
    % Create knots vector using bidirectional-average-accumulated-chord-length-parameterization method of
    % each isopar-line. Default knot interval is [0,1].

    % Input:
    %   A: Phy-coords of control points. size(A)=[3,numpoints].
    % Output:
    %   uv: Knot vector without 0 and 1 element. eg.[.1 .4 .5 .9].
    b=A(:,2:end)-A(:,1:end-1);
    db=sqrt(sum(b.*b));
    db=cumsum(db);% Cumulate sum.
    uv=db/sum(db);
    uv(end)=[]; % Remove the last element 1.
    end

end
%%
% srf=nrbtestsrf;
% bzrsrf=nrb2bzr(srf);
% bzr=bzrsrf{3,2};
% 
% [uv,lg]=flatcheck(bzr,.001);
% % [uv,lg]=flatcheck(bzr,10);
% 
% bzr.number=bzr.order;
% bzr.form='B-NURBS';
% 
% if lg==1
%     if iscell(uv)
%         for i=1:2
%             if ~isempty(uv{i})
%                 num=length(uv{i});
%                 for j=1:num
%                     crv{i,j}=nrbsrfextract(bzr, uv{i}(j), i);
%                 end
%             end
%         end
%     else
%         p=nrbeval(bzr,uv);
%     end
% end
% figure; nrbplot(srf,[50,50]);
% hold on;
% nrbplot(bzr,[50,50],'light','on','colormap','summer');
% numcrv=length(crv);
% for k=1:numcrv
%     nrbplot(crv{k}, 20);
% end
%                     



