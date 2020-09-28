function crvs = nrlextract(srf)

% NRBEXTRACT: construct NURL curves by extracting the boundaries 
%     of a NURL surface, or NURL surfaces by extracting the 
%     boundary of a NURL volume.
% 
% Calling Sequence:
% 
%   crvs = nrlextract(surf);
% 
% INPUT:
% 
%   surf        : NURL surface or volume, see nrbmak.
% 
% OUTPUT: 
% 
%   crvs        : array of NURL curves or NURL surfaces extracted.
% 
% Description:
% 
%  Constructs either an array of four NURL curves, by extracting the boundaries
%  of a NURL surface, or an array of six surfaces, by extracting the boundaries
%  of a NURL volume. The new entities are ordered in the following way
%
%    1: U = 0
%    2: U = 1
%    3: V = 0
%    4: V = 1
%    5: W = 0 (only for volumes)
%    6: W = 1 (only for volumes)
% 

if (~iscell (srf.knots))
  error('The boundary information is only extracted for NURL surfaces or volumes');  
end

if (numel (srf.knots) == 2)
  for ind = 1:2
    ind2 = mod (ind, 2) + 1;    % ind2 = [2 1];
    bnd1 = (ind - 1) * 2 + 1;
    bnd2 = (ind - 1) * 2 + 2;
    if (ind == 1)
      coefs1 = squeeze (srf.coefs(:,1,:));
      coefs2 = squeeze (srf.coefs(:,end,:));
    elseif (ind == 2)
      coefs1 = squeeze (srf.coefs(:,:,1));
      coefs2 = squeeze (srf.coefs(:,:,end));
    end
    crvs(bnd1) = nrlmake (coefs1, srf.knots{ind2}, srf.intervals{ind2}, srf.order(ind2));
    crvs(bnd2) = nrlmake (coefs2, srf.knots{ind2}, srf.intervals{ind2}, srf.order(ind2));
  end
elseif (numel (srf.knots) == 3)
  for ind = 1:3
    inds = setdiff (1:3, ind);
    bnd1 = (ind - 1) * 2 + 1;
    bnd2 = (ind - 1) * 2 + 2;
    if (ind == 1)
      coefs1 = squeeze (srf.coefs(:,1,:,:));
      coefs2 = squeeze (srf.coefs(:,end,:,:));
    elseif (ind == 2)
      coefs1 = squeeze (srf.coefs(:,:,1,:));
      coefs2 = squeeze (srf.coefs(:,:,end,:));
    elseif (ind == 3)
      coefs1 = squeeze (srf.coefs(:,:,:,1));
      coefs2 = squeeze (srf.coefs(:,:,:,end));
    end
    crvs(bnd1) = nrlmake (coefs1, {srf.knots{inds(1)} srf.knots{inds(2)}}, ...
        {srf.intervals{inds(1)} srf.intervals{inds(2)}}, [srf.order(inds(1)), srf.order(inds(2))]);
    crvs(bnd2) = nrlmake (coefs2, {srf.knots{inds(1)} srf.knots{inds(2)}}, ...
        {srf.intervals{inds(1)} srf.intervals{inds(2)}}, [srf.order(inds(1)), srf.order(inds(2))]);
  end
else
  error ('The entity is not a surface nor a volume')
end

end
