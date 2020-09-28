function [dist,npt,distava,nptava]=propt2line(pt,linept1,linept2)
% Get the distance from a given point and a line(defined by two vertices)
% and the coordinates of the projective point. This function is available
% for both 2D and 3D.

% Input:
%     pt: The given point.
%     linept1,linept2: Two vertices of the given line.
% Output:
%     dist: Distance from pt to line. t may not be located in available
%       domain [0,1]
%     npt: The projective point of pt. t may not be located in available
%       domain [0,1]
%     distava: Distance that is located in the available domain [0,1];
%     nptava: projectiive point that is located in [0,1].

% Note that the coordinate of the projection point is referred to the
% global coordinate system instead of the local coordinate.

slop=linept2-linept1;
% Represent the line using parametric representation. And solve the
% symbolic equation by using solve.m.
    syms t;
    lin=(1-t)*linept1+t*linept2;
    f=dot((lin-pt),slop);
    nt=solve(f==0);
    nt=eval(nt);
    % Don't consider the available domain, which means t can t<0 or t>1
    npt=subs(lin,t,nt);
    npt=eval(npt);
    % Constrain t in available domain, which means 0<=t<=1
    if t>=0 && t<=1
        nptava=npt;
    else
        if t<0
            nt=0;
        else
            nt=1;
        end
        nptava=subs(lin,t,nt);
        nptava=eval(nptava);
    end
  
% The projection point may locates out of the demain of the line.
% if length(slop)==2
%     if slop(1)~=0
%         k=slop(2)/slop(1);
%         b=linept1(2)-k*linept1(1);
%         npt(1,1)=(k*(pt(2)-b)+pt(1))/(k*k+1);
%         npt(1,2)=k*npt(1,1)+b;
%     else % The line has a infinite slope.
%         npt(1,1)=linept1(1);
%         npt(1,2)=pt(2);
%     end
% elseif length(slop)==3

dist=norm(npt-pt);
distava=norm(nptava-pt);
end

