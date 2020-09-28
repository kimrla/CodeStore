function [S, T, C]=getstc(tt, cc, di)

% Get the matrices of u and v coordinates for the triangle
% 
% Calling Sequences:
%
%     [S, T, C]=getstc(tt, cc, di)
%
% INPUTS:
% 
%      tt     - parametric evaluation points. The nurls should be  
%            a surface and tt is a cell {tu, tv} of the parametric coordinates
%
%     cc  - integration weights in one dimensional for u and v direction 
%             respectively, cc is a cell {cu, cv}.
%
%     di -  direction of integration 1 (u) or 2 (v)
%
% OUTPUT:
% 
%     S    :   a matrix of u coordinates for integration on a triangle
%     T    :   a matrix of v coordinates for integration on a triangle
%     C    :   a matrix of integration weights for triangle
% 

if iscell(tt)
    [v, u]=meshgrid(tt{2}, tt{1});
    [Cv, Cu]=meshgrid(cc{2}, cc{1}); 
    switch di
        case 1        
            S=u.*(1-v);
            T=v;       
            C=Cu.*Cv.*(1-v);
        case 2    
            S=u;
            T=v.*(1-u);
            C=Cu.*Cv.*(1-u);
    end
else
    error('The input parametric points should be a cell array.');
end





