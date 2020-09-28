function [nptslist1,nptslist2,npntslist]=intertrace(srf1,srf2,ptslist1,ptslist2,pntslist,varargin)
% Use tracing method to calculate the intersection points or correct the
% topology of intersection curve. Improved tracing method: if there is
% direction lost, use interbidi.m firstly to see if the 2 parts can be
% connected actually.

% Input:
%   srf1,srf2: NURBS/Bezier surface structure.
%   ptslist1,ptslist2: Par-coords of the PRECISE intersection ponit list of 
%       the 2 NURBS surfaces, which is a cell array, each cell represents 
%       a intersection part. numel(ptslist)=1 or 2.
%   pntslist: Phy-coords of the PRECISE intersection point list, 
%       corresponding to the ptslist1 or ptslist2.
%   varargin: Other input variables corresponding to 3/4-par iterative
%       (ite3par.m or ite4par.m), and stepvector.m.

% Output:
%   nptslist1,nptslist2: Par-coords of the PRECISE intersection ponit list of
%       the 2 NURBS surfaces AFTER using the tracing method to update.
%   npntslist: Phy-coords corresponding to ptslist1 or ptslist2.

% lg=0: Doesn't arrive to the boundary or one specific point, needing to
%   continue to trace, or transform to use ite3par.
% lg=1: Has arrived to the boudnary or the point or 'enter' the intersection
%   list (meaning the intersection curve is closed or connected). Stop to trace.
% lg=-1: Distance d doesn't satisfy the precision, meaning the intersection
%   point is NOT calculated normally.

% Three tracing ending condition determination:
% 1.traceend_distance.m
% 2.traceend_boundary.m
% 3.traceend_point.m

%注意1：被调函数中如果有可变输入，要把所有的输入写全
%注意2：传递函数句柄时，要注意是否真的是引用传递，而不是每次只传递了一个确定的数

% Default value of the optional variance.
ITE=5;
MIU=[];
TOL=1e-10;
STEP=1E-2;
ANGLE=2*pi/1e4;
MODE='fix';
TOLITE34=0.1; % Tolerance of transforming from ite4par.m to ite3par.m.

options=struct('iterative',ITE,...
    'mountain',MIU,...
    'tolerance',TOL,...
    'step',STEP,...
    'angle',ANGLE,...
    'mode',MODE,...
    'tolerance34',TOLITE34);
optionnames=fieldnames(options);

nargin=length(varargin);
if round(nargin/2)~=nargin/2
    error('The input is error!');
end

for pair=reshape(varargin,2,[])
    if any(strcmp(pair{1},optionnames))
        options.(pair{1})=pair{2};
    end
end

nump=numel(ptslist1);
if nump==1% Input only one intersection part.
          
    for ii=1:2% Check the 2 ending points respectively
            [ptslist1,ptslist2,pntslist,pts_end1,pts_end2,pts_end,...
                pnts_end]=ptsdeal(ii,ptslist1,ptslist2,pntslist);


            % Preliminarily check whether the points are located in or near the
            % boundary. If they are ABSOLUTELY NOT, then aboid to call the
            % corresponding functions traceend_.m.
            
            [~,~,neart]=predistance(pts_end([2,4],:),options.step);
            while true
                if isempty(neart) % Points are located in or near the boundary
                
                    [bdlg,bduv,bdst]=traceend_boundary(pts_end1(2),...
                            pts_end2(2),options.tolerance,options.tolerance34);

                    if bdlg==1% Ending point is located in the boundary
                        break;
                    else% Ending point is NOT located in but near the boundary.
                        if ~isempty([bduv{:}])% Select srf1 as the iso-par-line
                            if ~isempty(bduv{1})
                                signal=1;% Glue the 1st par(u).
                                signal_=2;
                                st2=bduv{1};
                            else
                                signal=2;% Glue the 2nd par(v).
                                signal_=1;
                                st2=bduv{1};
                            end
                            st2(1)=[];
                            % 3 par iterative ite3par.m
                            for jj=1:numel(st2)
                                sv=stepvector(srf1,srf2,pts_end1{2},...
                                    pts_end2{2},pnts_end{2},...
                                    'step',abs(st2(2)-st2(1)),...
                                    'angle',options.angle,...
                                    'tolerance',options.tolerance,...
                                    'mode',options.mode);
                                % Need to check if the direction of the step
                                % vector is correct or oppposite.
                                newpnts0=pnts_end{2}+sv;
                                if dot(newpnts0-pntslist(end,:),...
                                        pntslist(end-1,:)-pntslist(end,:))>0
                                    sv=-sv;
                                    newpnts0=pnts_end{2}+sv;
                                end

                                [nextpts,nextpnts,nextd]=ite3par(srf2,srf1,...,
                                    [signal,st2(jj)],[pts_end2{2},...
                                    pts_end1{2}(signal_)],newpnts0,...
                                    'iterative',options.iterative,...
                                    'mountain',options.mountain,...
                                    'tolerance',options.tolerance);
                                if traceend_distance(nextd,options.tolerance)==-1
                                    % The new intersection point is NOT
                                    % obtained successfully.                             
                                    break;% Stop this ending point tracing   
                                else
                                    % Update the intersection list and continue
                                    % the circulation.
                                    ptslist1(end+1,:)=nextpts(1,:);
                                    ptslist2(end+1,:)=nextpts(2,:);
                                    pntslist(end+1,:)=(nextpnts(1,:)+...
                                        nextpnts(2,:))/2;
                                    [ptslist1,ptslist2,pntslist,pts_end1,...
                                        pts_end2,pts_end,pnts_end]=ptsdeal(ii,...
                                        ptslist1,ptslist2,pntslist);                
                                end  
                            end                           
                            break;% Break the circulation while true.
                                                      
                        elseif ~isempty([bdst{:}])% Select srf2 as the iso-par-line
                            if ~isempty(bdst{1})
                                signal=1;% Glue the 1st par(s).
                                st2=bdst{1};
                            else
                                signal=1;% Glue the 2nd par(t).
                                st2=bdst{1};
                            end
                            st2(1)=[];
                             % 3 par iterative ite3par.m
                            for jj=1:numel(st2)
                                sv=stepvector(srf1,srf2,pts_end1{2},...
                                    pts_end2{2},pnts_end{2},...
                                    'step',abs(st2(2)-st2(1)),...
                                    'angle',options.angle,...
                                    'tolerance',options.tolerance,...
                                    'mode',options.mode);
                                % Need to check if the direction of the step
                                % vector is correct or oppposite.
                                newpnts0=pnts_end{2}+sv;
                                if dot(newpnts0-pntslist(end,:),...
                                        pntslist(end-1,:)-pntslist(end,:))>0
                                    sv=-sv;
                                    newpnts0=pnts_end{2}+sv;
                                end

                                [nextpts,nextpnts,nextd]=ite3par(srf1,srf2,...,
                                    [signal,st2(jj)],[pts_end1{2},...
                                    pts_end2{2}(signal_)],newpnts0,...
                                    'iterative',options.iterative,...
                                    'mountain',options.mountain,...
                                    'tolerance',options.tolerance);
                                if traceend_distance(nextd,options.tolerance)==-1
                                    % The new intersection point is NOT
                                    % obtained successfully.                             
                                    break;% Stop this ending point tracing   
                                else
                                    % Update the intersection list and continue
                                    % the circulation.
                                    ptslist1(end+1,:)=nextpts(1,:);
                                    ptslist2(end+1,:)=nextpts(2,:);
                                    pntslist(end+1,:)=(nextpnts(1,:)+...
                                        nextpnts(2,:))/2;
                                    [ptslist1,ptslist2,pntslist,pts_end1,...
                                        pts_end2,pts_end,pnts_end]=ptsdeal(ii,...
                                        ptslist1,ptslist2,pntslist);                
                                end  
                            end
                            
                            break;% Break the circulation while true.

                        else % Near the boundary but use ite4par.m
                            sv=stepvector(srf1,srf2,pts_end1{2},...
                                pts_end2{2},pnts_end{2},...
                                'step',options.step,...
                                'angle',options.angle,...
                                'tolerance',options.tolerance,...
                                'mode',options.mode);
                            % Need to check if the direction of the step
                            % vector is correct or oppposite.
                            newpnts0=pnts_end{2}+sv;
                            if dot(newpnts0-pntslist(end,:),...
                                    pntslist(end-1,:)-pntslist(end,:))>0
                                sv=-sv;
                                newpnts0=pnts_end{2}+sv;
                            end

                            [nextpts,nextpnts,nextd]=ite4par(srf1,srf2,...,
                                [pts_end1{2},pts_end2{2}],newpnts0,...
                                'iterative',options.iterative,...
                                'mountain',options.mountain,...
                                'tolerance',options.tolerance);
                            if traceend_distance(nextd,options.tolerance)==-1
                                % The new intersection point is NOT
                                % obtained successfully.                             
                                break;% Break the circulation while true.  
                            else
                                % Update the intersection list and continue
                                % the circulation.
                                ptslist1(end+1,:)=nextpts(1,:);
                                ptslist2(end+1,:)=nextpts(2,:);
                                pntslist(end+1,:)=(nextpnts(1,:)+...
                                    nextpnts(2,:))/2;
                                [ptslist1,ptslist2,pntslist,pts_end1,...
                                    pts_end2,pts_end,pnts_end]=ptsdeal(ii,...
                                    ptslist1,ptslist2,pntslist);                
                            end  
                        end
                    end
                    % If use ite3par, then break the circulation. If use
                    % ite4par, then continue the circulation.
                    if ~isempty([bduv{:}]) || ~isempty([bdst{:}])
                        break;
                    end
                                      
                else% Points are far away from the boundary. Check point traceend_point.
                    % Check if the curve has extended to the 1st ending point.
                    [~,neard1,neart1]=predistance(pts_end1{1},options.step,pts_end1{2});
                    [~,neard2,neart2]=predistance(pts_end2{1},options.step,pts_end2{2});
                    [neard,nearid]=min(neard1,neard2);

                    if isempty(neart1) || isempty(neart2)% Ending point is located in or near the 1st point.
                        if nearid==1% Deal with srf1
                            [lg,uv]=traceend_point(pts_end1{2},pnts_end{2},...
                                ptslist1,pntslist,...
                                options.tolerance,options.tolerance34);
                            if lg==1% Point is located in the ending point.
                                % If the curve is closed, the 1st element is
                                % the same as the last one( in precision error)
                                nptslist1{1}=ptslist1;
                                nptslist2{1}=ptslist2;
                                npntslist{1}=npntslist;
                                return;
                            else% Points are near the ending point, use ite3par.
                                if abs(uv{1}(end)-uv{1}(1))>abs(uv{2}(end)-uv{2}(1))
                                    st2=uv{2};
                                    st2(1)=[];
                                    signal=2;
                                    signal_=1;
                                else
                                    st2=uv{1};
                                    st2(1)=[];
                                    signal=1;
                                    signal_=2;
                                end

                                for jj=1:numel(st2)
                                    sv=stepvector(srf1,srf2,pts_end1{2},...
                                        pts_end2{2},pnts_end{2},...
                                        'step',abs(st2(2)-st2(1)),...
                                        'angle',options.angle,...
                                        'tolerance',options.tolerance,...
                                        'mode',options.mode);
                                    % Need to check if the direction of the step
                                    % vector is correct or oppposite.
                                    newpnts0=pnts_end{2}+sv;
                                    if dot(newpnts0-pntslist(end,:),...
                                            pntslist(end-1,:)-pntslist(end,:))>0
                                        sv=-sv;
                                        newpnts0=pnts_end{2}+sv;
                                    end

                                    [nextpts,nextpnts,nextd]=ite3par(srf2,srf1,...,
                                        [signal,st2(jj)],[pts_end2{2},...
                                        pts_end1{2}(signal_)],newpnts0,...
                                        'iterative',options.iterative,...
                                        'mountain',options.mountain,...
                                        'tolerance',options.tolerance);
                                    if traceend_distance(nextd,options.tolerance)==-1
                                        % The new intersection point is NOT
                                        % obtained successfully.                             
                                        break;% Stop this ending point tracing   
                                    else
                                        % Update the intersection list and continue
                                        % the circulation.
                                        ptslist1(end+1,:)=nextpts(1,:);
                                        ptslist2(end+1,:)=nextpts(2,:);
                                        pntslist(end+1,:)=(nextpnts(1,:)+...
                                            nextpnts(2,:))/2;
                                        [ptslist1,ptslist2,pntslist,pts_end1,...
                                            pts_end2,pts_end,pnts_end]=ptsdeal(ii,...
                                            ptslist1,ptslist2,pntslist);                
                                    end  
                                end  
                                break;% Break the circulation while true.
                            end

                        else% Deal with srf2
                            [lg,uv]=traceend_point(pts_end2{2},pnts_end{2},...
                                ptslist2,pntslist,...
                                options.tolerance,options.tolerance34);   
                            if lg==1% Point is located in the ending point.
                                % If the curve is closed, the 1st element is
                                % the same as the last one( in precision error)
                                nptslist1{1}=ptslist1;
                                nptslist2{1}=ptslist2;
                                npntslist{1}=npntslist;
                                return;
                            else% Points are near the ending point, use ite3par.
                                if abs(uv{1}(end)-uv{1}(1))>abs(uv{2}(end)-uv{2}(1))
                                    st2=uv{2};
                                    st2(1)=[];
                                    signal=2;
                                    signal_=1;
                                else
                                    st2=uv{1};
                                    st2(1)=[];
                                    signal=1;
                                    signal_=2;
                                end

                                for jj=1:numel(st2)
                                    sv=stepvector(srf1,srf2,pts_end2{2},...
                                        pts_end1{2},pnts_end{2},...
                                        'step',abs(st2(2)-st2(1)),...
                                        'angle',options.angle,...
                                        'tolerance',options.tolerance,...
                                        'mode',options.mode);
                                    % Need to check if the direction of the step
                                    % vector is correct or oppposite.
                                    newpnts0=pnts_end{2}+sv;
                                    if dot(newpnts0-pntslist(end,:),...
                                            pntslist(end-1,:)-pntslist(end,:))>0
                                        sv=-sv;
                                        newpnts0=pnts_end{2}+sv;
                                    end

                                    [nextpts,nextpnts,nextd]=ite3par(srf1,srf2,...,
                                        [signal,st2(jj)],[pts_end1{2},...
                                        pts_end2{2}(signal_)],newpnts0,...
                                        'iterative',options.iterative,...
                                        'mountain',options.mountain,...
                                        'tolerance',options.tolerance);
                                    if traceend_distance(nextd,options.tolerance)==-1
                                        % The new intersection point is NOT
                                        % obtained successfully.                             
                                        break;% Stop this ending point tracing   
                                    else
                                        % Update the intersection list and continue
                                        % the circulation.
                                        ptslist1(end+1,:)=nextpts(1,:);
                                        ptslist2(end+1,:)=nextpts(2,:);
                                        pntslist(end+1,:)=(nextpnts(1,:)+...
                                            nextpnts(2,:))/2;
                                        [ptslist1,ptslist2,pntslist,pts_end1,...
                                            pts_end2,pts_end,pnts_end]=ptsdeal(ii,...
                                            ptslist1,ptslist2,pntslist);                
                                    end  
                                end  
                                break;% Break the circulation while true.
                            end   
                        end
                        
                    else% Both srf1 and srf2 are ordinary ending points.Use 
                        % ite4par and consider neart.
                        if nearid==1
                            neart=neart1;                      
                        else
                            neart=neart2;                       
                        end
                        for kk=1:neart
                            sv=stepvector(srf1,srf2,pts_end1{2},...
                                pts_end2{2},pnts_end{2},...
                                'step',options.step,...
                                'angle',options.angle,...
                                'tolerance',options.tolerance,...
                                'mode',options.mode);
                            % Need to check if the direction of the step
                            % vector is correct or oppposite.
                            newpnts0=pnts_end{2}+sv;
                            if dot(newpnts0-pntslist(end,:),...
                                    pntslist(end-1,:)-pntslist(end,:))>0
                                sv=-sv;
                                newpnts0=pnts_end{2}+sv;
                            end

                            [nextpts,nextpnts,nextd]=ite4par(srf1,srf2,...,
                                [pts_end1{2},pts_end2{2}],newpnts0,...
                                'iterative',options.iterative,...
                                'mountain',options.mountain,...
                                'tolerance',options.tolerance);
                            if traceend_distance(nextd,options.tolerance)==-1

                                break;% Stop this ending point tracing and 
                                % check another ending point
                            else
                                % Update the intersection list and continue
                                % the circulation.
                                ptslist1(end+1,:)=nextpts(1,:);
                                ptslist2(end+1,:)=nextpts(2,:);
                                pntslist(end+1,:)=(nextpnts(1,:)+...
                                    nextpnts(2,:))/2;
                                [ptslist1,ptslist2,pntslist,pts_end1,...
                                    pts_end2,pts_end,pnts_end]=ptsdeal(ii,...
                                    ptslist1,ptslist2,pntslist);                
                            end  
                        end
                        if kk==neart% Return the traceend_boundary
                            continue;
                        else% Return for ii-1:2 to check another ending point.
                            break;
                        end
                    end               
                end
            end
        end
      
    nptslist1{1}=ptslist1;
    nptslist2{1}=ptslist2;
    npntslist{1}=npntslist;

    
elseif nump==2 % Input 2 intersection parts.
    
    
    
else
    error('The input intersection points should be 1 or 2 cell array!');
end





end



function [id,d,t]=predistance(pts,step,varargin)
% Preliminarily check whether the point is far away from the 4 boundaries
% in par-domains of the 2 surfaces or one specific point 
% (testpt=varargin{1}) . If the point is far, then export the
% nearest point INDEX and its distance (d) to the boundary and the times(t) 
% that does NOT need to call the functions traceend_.m.
% size(pts)=[nump,2]. step is the step length of trace method.
% If the points are NOT far away from the boundaries, then t=[].
if nargin==3
    testpt=varargin{1};
else
    testpt=[0.1,0.9];% Default boundary value
end
N=10; % Tolerance
n=size(pts,1);
d=[abs(pts(:)-testpt(1));abs(pts(:)-testpt(2))];
[d,id]=min(d);
id=rem(id,n);
if id==0
    id=n;
end
t=fix(d/step);
if t<N
    % If the distance exceeeds more than N times, regard this point as
    % ordinary point. The distance is from point to 0.1 or 0.9.
    t=[];
end

end

function [ptslist1,ptslist2,pntslist,pts_end1,pts_end2,pts_end,pnts_end]=ptsdeal(ii,ptslist1,ptslist2,pntslist)
if iscell(ptslist1)
    ptslist1=ptslist1{1};
    ptslist2=ptslist2{1}; 
    pntslist=pntslist{1}; 
    % Make sure size(pts)=[numpoints,2]
    ptslist1=shiftdim(ptslist1,find(size(ptslist1)==2));
    % Make sure size(pts)=[numpoints,2]
    ptslist2=shiftdim(ptslist2,find(size(ptslist2)==2));
    % Make sure size(pts)=[numpoints,3]
    pntslist=shiftdim(pntslist,find(size(pntslist)==3));
end
if ii==2% Reverse the direction of the intersection list
    ptslist1=ptslist1(end:1,:);
    ptslist2=ptslist2(end:1,:);
    pntslist=pntslist(end:1,:);
end

% 2 ending vertices of the input intersection curve part.
pts_end1{1}=ptslist1(1,:);
pts_end1{2}=ptslist1(end,:); % On srf1
pts_end2{1}=ptslist2(1,:);
pts_end2{2}=ptslist2(end,:); % On srf2
pnts_end{1}=pntslist(1,:);% 1st ending point in the list
pnts_end{2}=pntslist(end,:);% last ending point in the list
% size(pts_end)=[4,2]
pts_end=[pts_end1{1};pts_end1{2};pts_end2{1};pts_end2{2}];
% 'Extend' the curve part always from the last ending point {2}. 



end






