function tsrf = nrltransp(srf)
% 
% NRLTRANSP: Transpose a NURL surface, by swapping U and V directions.
% 
% Calling Sequence:
% 
%   tsrf = nrltransp(srf)
%
% INPUT:
% 
%   srf		: NURL surface, see nrbmak.
%
% OUTPUT:
% 
%   tsrf	: NURL surface with U and V diretions transposed.
% 
% Description:
% 
%   Utility function that transposes a NURL surface, by swapping U and
%   V directions. NURL curves cannot be transposed.
%

if ~iscell(srf.knots)
  error(' A NURL curve cannot be transposed.');
elseif size(srf.knots,2) == 3
  error('The transposition of NURL volumes has not been implemented.');
end  

tsrf = nrlmake(permute(srf.coefs,[1 3 2]), fliplr(srf.knots), fliplr(srf.intervals), fliplr(srf.order));

end

%% demo - 1
%  srf = nrl4surf([0 0 0], [1 0 1], [0 1 1], [1 1 2]);
%  nrlplot(srf,[20 5]);
%  title('Plane surface and its transposed (translated)')
%  hold on
%  srf.coefs(3,:,:) = srf.coefs(3,:,:) + 5;
%  srf = nrltransp(srf);
%  nrlplot(srf,[20 5]);
%  hold off

%% demo - 2
%  srf  = nrlrevolve(nrlline([1 0],[2 0]), [0 0 0], [0 0 1], pi/2);
%  srft = nrltransp(srf);
%  nrlplot(srf,[20 5]);
%  figure; 
%  nrlplot(srft,[5 20]);





