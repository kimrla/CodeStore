function crv=nrlsrfextract(srf, s, t)

% nrlcrvextract: extract a curve from a nurl surface
%
% Calling Sequence:
% 
%    crv=nrlsrfextract(srf, [], t)
%
%    crv=nrlsrfextract(srf, s, [])
% 
% INPUT:
% 
%    srf   :  Anurl surfaces, see nrlmake.
%
%    s, t   :   Pamatric point. One of them should be empty.
%               Curve with respect to the nonempty one will be extracted.
% 
% OUTPUT:
%
%     crv - a nurl curve extracted from the surface
%

if isempty(s) && ~isempty(t)
    [pts, w]=nrleval(srf, [srf.knots{1}; repmat(t, 1, srf.number(1))]);
    pts=pts./repmat(w, 3, 1);
    coefs=[pts; w];
    crv=nrlmake(coefs, srf.knots{1}, srf.intervals{1}, srf.order(1));
elseif ~isempty(s) && isempty(t)
    [pts, w]=nrleval(srf, [repmat(s, 1, srf.number(2)); srf.knots{2}]);
    pts=pts./repmat(w, 3, 1);
    coefs=[pts; w];
    crv=nrlmake(coefs, srf.knots{2}, srf.intervals{2}, srf.order(2));
end

%% demo
% % Create a nurl surfaces
% srf = nrltestsrf; 
% crv1=nrlsrfextract(srf, [], 0.2);
% crv2=nrlsrfextract(srf, 0.2, []);
% 
% % Plot the surface
% figure; hold on; 
% nrlplot(srf, [100, 100]); 
% nrlplot(crv1, 100); 
% nrlplot(crv2, 100); 
% axis equal; view(3); 
% shading interp; 
