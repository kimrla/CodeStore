function bzr=nrb2bzr(nrb,varargin)

% Transform the NURBS curve or surface to Bezier model. 

% Input:
%   nrb: NURBS model. If the NURBS model is exactly the Bezier model, then
%       the output bzr contains the original information of the input.
%   uv: Par-cord of the decomposed point. If there's no uv, the the default 
%       decomposition is implemented in all the break points.If there are 
%       several new points, then uv=[u1,u2...] for a curve's par-domain; uv={[u1,u2...um],[v1,v2...vn]} for a surface domain.
% Output:
%   bzr: Cell array of a Bezier model structure, which stores all the Bezier segments.

% The Bezier structure: 
%   Bezier.form: 'Bezier'
%   Bezier.dim: Dimension of the control points. dim=3(inhomogeneous) or dim=4(homogeneous).
%   Bezier.order: Order of the model.
%   Bezier.coefs: Control points' cords. If it's the rational Bezier, then
%       the coefs contains homogeneous cords. size(coefs)=[4,numpoints] or [3,numpoints].
%   Bezier.knots: Local cords of the original knots. The transform from
%       local cords to global cords [0,1] is t=(1-u)*a+u*b, where
%       t=[a,b],u=[0,1]. The repeatability of each interval equals order,
%       which means all Bezier segments are clamped.
%   Bezier.nrb: The original NURBS structure before the transformation.

% The new several surfaces are numbered as vector form:
% bzr{1,1} bzr{1,2} bzr{1,3}
% bzr{2,1} bzr{2,2} bzr{2,3}
% bzr{3,1} bzr{3,2} bzr{3,3}

deg=nrb.order-1;
dim=length(deg);
knt=nrb.knots;
newuv=[];
if dim==2 
    newuv={[],[]} ; 
end

if nargin>1
    newuv=varargin{1};
end
if dim==1% Bezier curve
    brk=unique([knt,newuv]);
    newknt=kntbrkdegmult(brk,deg,deg);
    kntins=new_knots(knt,newknt);% Find points needed to be insterted with their repeatability.
    newnrb=nrbkntins(nrb,kntins);
    % Extract each Bezier segment after the transformation. Add the field bzr.coefs and bzr.knots.
    numbrk=length(brk);
    pts=newnrb.coefs;
    j=1;
    for i=1:numbrk-1
        knots{i}=[repmat(brk(i),1,deg+1),repmat(brk(i+1),1,deg+1)];
        coefs{i}=pts(:,j:j+deg);
        j=j+deg;
        bzr{i}.knots=knots{i};
        bzr{i}.coefs=coefs{i};
        bzr{i}.form='Bezier';
        bzr{i}.nrb=nrb;
        bzr{i}.dim=size(pts,1);
        bzr{i}.order=newnrb.order;
%         bzr{i}.number= bzr{i}.order
    end
%     bzr{1}.knots=[0,knots{1}];bzr{end}.knots=[knots{end},1];
elseif dim==2% Bezier surface
    for i=1:2
        brk{i}=unique([knt{i},newuv{i}]);
    end
    % Use the cellfun to reduce the code
    newknt=cellfun(@kntbrkdegmult,brk,num2cell(deg),num2cell(deg),'UniformOutput',false);
    kntins=cellfun(@new_knots,knt,newknt,'UniformOutput',false);
    newnrb=nrbkntins(nrb,kntins);
    % Extract each Bezier segment after the transformation. Add the field bzr.coefs and bzr.knots.
    numbrk(1)=length(brk{1});numbrk(2)=length(brk{2});
    pts=newnrb.coefs;
    f=1;g=1;
    for j=1:numbrk(2)-1  % Number of rows: number of par-points in v-direction
        for i=1:numbrk(1)-1 % Number of columns: number of par-points in u-direction
            knots{j,i}={[repmat(brk{1}(i),1,deg(1)+1),repmat(brk{1}(i+1),1,deg(1)+1)],[repmat(brk{2}(j),1,deg(2)+1),repmat(brk{2}(j+1),1,deg(2)+1)]};
            coefs{j,i}=pts(:,f:f+deg(1),g:g+deg(2));
            f=f+deg(1);
            bzr{j,i}.knots=knots{j,i};
            bzr{j,i}.coefs=coefs{j,i};
            bzr{j,i}.form='Bezier';
            bzr{j,i}.nrb=nrb;
            bzr{j,i}.dim=size(pts,1);
            bzr{j,i}.order=newnrb.order;
%             bzr{j,i}.number=bzr{j,i}.order;
        end
        g=g+deg(2);
        f=1;
    end
%     for j=1:numbrk(2)-1
%         bzr{j,1}.knots{1}=[0,knots{j,1}{1}];bzr{j,end}.knots{1}=[knots{j,end}{1},1];
%     end
%     for i=1:numbrk(1)-1
%         bzr{1,i}.knots{2}=[0,knots{1,i}{2}];bzr{end,i}.knots{2}=[knots{end,i}{2},1];
%     end

end

end

%% demo
% % Need to create the plot function to draw the BEZIER model.
% % bzr.form~='NURBS'
% % There is no bzr.number, which equals bzr.order.
% % Each knot interval's repeatability is error because each segment is not clamped.

% pnts = zeros(3,5,5);
% pnts(:,:,1) = [ 0.0  3.0  5.0  8.0 10.0;     % w*x
%                 0.0  0.0  0.0  0.0  0.0;     % w*y
%                 2.0  2.0  7.0  7.0  8.0];    % w*z
% 
% pnts(:,:,2) = [ 0.0  3.0  5.0  8.0 10.0;
%                 3.0  3.0  3.0  3.0  3.0;
%                 0.0  0.0  5.0  5.0  7.0];
% 
% pnts(:,:,3) = [ 0.0  3.0  5.0  8.0 10.0;
%                 5.0  5.0  5.0  5.0  5.0;
%                 0.0  0.0  5.0  5.0  7.0];
% 
% pnts(:,:,4) = [ 0.0  3.0  5.0  8.0 10.0;
%                 8.0  8.0  8.0  8.0  8.0;
%                 5.0  5.0  8.0  8.0 10.0];
% 
% pnts(:,:,5) = [ 0.0  3.0  5.0  8.0 10.0;
%                10.0 10.0 10.0 10.0 10.0;
%                 5.0  5.0  8.0  8.0 10.0];
%             
% pnts(:,4,:)=[];
% 
% % knots
% knots{1} = [0 0 0 1/2 1 1 1]; % knots along u
% knots{2} = [0 0 0 1/3 2/3 1 1 1]; % knots along v
% 
% % make and draw nurbs surface
% srf = nrbmak(pnts,knots);
% figure;
% hold on;
% % nrbctrlplot(srf);
% nrbplot(srf,[50,50]);
% axis equal;
% 
% % [~,numcol,numrow]=size(pnts);
% % for i=1:numrow
% %     for j=1:numcol
% %         text(pnts(1,j,i),pnts(2,j,i),pnts(3,j,i),{i,j});
% %     end
% % end
%    
% bzrsrf=nrb2bzr(srf);
% [m,n]=size(bzrsrf);
% bzr=bzrsrf{2,2};
% bzr.form='B-NURBS';bzr.number=bzr.order;
% nrbctrlplot(bzr);

