function nbzr=bzrseg(bzr,uv)

% Segment a Bezier model(curve or surface) using the input par-points. The
% input points can be more than one.

% Input:
%   bzr: Bezier structure. See nrb2bzr.m.
%   uv: Parameter coordinates of segment points. If uv is a cell array,
%       then uv={[u],[v]} stores u and v cords respectively of the Bezier
%       surface. (eg.uv={[0,0,0.4],[]}).
% Output:
%   nbzr: New Bezier models after segment, which is a cell array.

bzr.number=bzr.order;
nbzr=nrb2bzr(bzr,uv);

end
   

    
