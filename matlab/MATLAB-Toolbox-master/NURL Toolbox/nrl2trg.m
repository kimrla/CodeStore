function [tsrf, du]=nrl2trg(srf)

% Transform a nurl triangle patch for subsequent manipulations. 
% 
%  Discription : This function also check whether the pacth is 
%     an ellipse, quadrangle, or triangle 
% 
% Calling Sequences:
%
%     [tsrf, du]=nrl2trg(srf)
%
% INPUTS:
%
%      srf - a nurl surface
%
% OUTPUT:
% 
%     tsrf    :   a nurl triangle. If du<0, nothing is done
%                   and an empty variable is returned
%     du :   an index of the shape  of the patch or the direction of u axis
%               du = -1 means the patch is an ellipse 
%               du = 0 means the patch is a quadrangle 
%               du = 1 means the patch is a triangle and v axis
%                        is transformed to point to the focus point 
% 

% Transform the surface to make v axis point to focus point
[srf, du]=trgsrfdirect(srf); 
if du<0
    tsrf=[];
    return;
end

% Extract the three boundary curves
crvt = nrlextract(srf);
crvs=crvt([3, 2, 1]);

% Make a triangle nurls patch
tsrf.form='T-NURL';
tsrf.faces=srf;
tsrf.edges=crvs;




