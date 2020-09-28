function srf = nrlcoons(u1, u2, varargin)
% 
% NRLCOONS: Construction of a Coons patch.
% 
% Calling Sequence:
% 
%   srf = nrlcoons(ucrv1, ucrv2, vcrv1, vcrv2)
% 
% INPUT:
% 
%   ucrv1	: NURL curve defining the bottom U direction boundary of
% 		the constructed NURL surface.
% 
%   ucrv2	: NURL curve defining the top U direction boundary of
% 		the constructed NURL surface.
% 
%   vcrv1	: NURL curve defining the bottom V direction boundary of
% 		the constructed NURL surface.
% 
%   vcrv2	: NURL curve defining the top V direction boundary of
% 		the constructed NURL surface.
%
% OUTPUT:
% 
%   srf		: Coons NURL surface patch.
% 
% Description:
% 
%   Construction of a bilinearly blended Coons surface patch from four NURL
%   curves that define the boundary.
% 
%   The orientation of the four NURL boundary curves.
% 
%          ^ V direction
%          |
%          |     ucrv2
%          ------->--------
%          |              |
%          |              |
%    vcrv1 ^   Surface    ^ vcrv2
%          |              |
%          |              |
%          ------->-----------> U direction
%                ucrv1
%

narg=length(varargin);
if narg>2
    error('Incorrect number of input arguments'); 
elseif narg ==0
    srf = nrlruled (u1, u2);
elseif narg ==1
    v1=varargin{1};
    srf = nrltrgcoons(u1, u2, v1);
elseif narg ==2
    v1=varargin{1};
    v2=varargin{2};
    nrls={u1, u2, v1, v2};
    lth=[nrlmeasure(nrls{1})
       nrlmeasure(nrls{2})
       nrlmeasure(nrls{3})
       nrlmeasure(nrls{4}) ];
    ii=find(lth>(1e-3)*max(lth));
    if length(find(ii))==3
        nrls=nrls(ii);
        crv1=nrls{1}; 
        crv2=nrls{2}; 
        crv3=nrls{3}; 
        srf = nrltrgcoons(crv1, crv2, crv3);
    elseif length(find(ii))==2
        nrls=nrls(ii);
        crv1=nrls{1}; 
        crv2=nrls{2}; 
        srf = nrlruled (crv1, crv2);
    else
        tnrls=trans4coosnsrf(nrls);
        u1=tnrls{1}; u2=tnrls{2};
        v1=tnrls{3}; v2=tnrls{4}; 

        if (max (abs (nrleval (u1, u1.knots(1)) - nrleval (v1, v1.knots(1)))) > 1e-10 || ...
            max (abs (nrleval (u1, u1.knots(end)) - nrleval (v2, v2.knots(1)))) > 1e-10 || ...
            max (abs (nrleval (u2, u2.knots(1)) - nrleval (v1, v1.knots(end)))) > 1e-10 || ...
            max (abs (nrleval (u2, u2.knots(end)) - nrleval (v2, v2.knots(end)))) > 1e-10)
          error ('The four curves do not define a closed boundary')
        end

        r1 = nrlruled(u1, u2);
        r2 = nrltransp(nrlruled(v1, v2));
        t  = nrl4surf(u1.coefs(1:3,1), u1.coefs(1:3,end), u2.coefs(1:3,1), u2.coefs(1:3,end));

        % raise all surfaces to a common degree
        du = max([r1.order(1), r2.order(1), t.order(1)]);
        dv = max([r1.order(2), r2.order(2), t.order(2)]);
        r1 = nrldegelev(r1, [du - r1.order(1), dv - r1.order(2)]);
        r2 = nrldegelev(r2, [du - r2.order(1), dv - r2.order(2)]);
        t  = nrldegelev(t,  [du - t.order(1),  dv - t.order(2)]);

        % merge the knot vectors, to obtain a common knot vector

        % U intervals
        k1 = r1.intervals{1};
        k2 = r2.intervals{1};
        k3 = t.intervals{1};
        ku = unique([k1 k2 k3]);

        % V intervals
        k1 = r1.intervals{2};
        k2 = r2.intervals{2};
        k3 = t.intervals{2};
        kv = unique([k1 k2 k3]);

        % Insert intervals
        r1 = nrlintins(r1, {ku, kv});
        r2 = nrlintins(r2, {ku, kv});
        t  = nrlintins(t,  {ku, kv});

        % U knots
        nu1=countknts(r1.intervals{1}, r1.knots{1}); 
        nu2=countknts(r2.intervals{1}, r2.knots{1});
        nu3=countknts(t.intervals{1}, t.knots{1});
        nu=max([nu1; nu2; nu3]);
        inu1=nu-nu1; inu2=nu-nu2; inu3=nu-nu3;

        % V knots
        nv1=countknts(r1.intervals{2}, r1.knots{2}); 
        nv2=countknts(r2.intervals{2}, r2.knots{2});
        nv3=countknts(t.intervals{2}, t.knots{2});
        nv=max([nv1; nv2; nv3]);
        inv1=nv-nv1; inv2=nv-nv2; inv3=nv-nv3;

        % Insert knots
        r1 = nrlkntins(r1, {inu1, inv1});
        r2 = nrlkntins(r2, {inu2, inv2});
        t  = nrlkntins(t,  {inu3, inv3});

        % combine coefficient to construct Coons surface
        coefs(1,:,:) = r1.coefs(1,:,:) + r2.coefs(1,:,:) - t.coefs(1,:,:);
        coefs(2,:,:) = r1.coefs(2,:,:) + r2.coefs(2,:,:) - t.coefs(2,:,:);
        coefs(3,:,:) = r1.coefs(3,:,:) + r2.coefs(3,:,:) - t.coefs(3,:,:);
        coefs(4,:,:) = r1.coefs(4,:,:) + r2.coefs(4,:,:) - t.coefs(4,:,:);
        srf = nrlmake(coefs, r1.knots, r1.intervals, r1.order);
    end
end

end

% Transform 4 nurl curves for a coons surface
function tnrls=trans4coosnsrf(nrls)

    % Take nrls{1} as reference curve 1
    tnrls = cell(1,4);
    tnrls{1}=nrls{1};    
    [xyz1, r1] = nrlaxis(nrls{1});
    nrls(1) = [];
    toler=(5e-3)*norm(r1);

    % Find curve 3 that have the same origin as curve 1
    for j=1:3
        rnrlj=nrlreverse(nrls{j});
        xyz = nrlaxis(nrls{j});
        rxyz = nrlaxis(rnrlj);
        if norm(xyz-xyz1)<toler
            tnrls{3}=nrls{j};
            nrls(j) = []; 
            break;
        elseif norm(rxyz-xyz1)<toler
            tnrls{3}=rnrlj;
            nrls(j) = []; 
            break;
        end
    end

    % Find curve 4 whose origin is the end of curve 1
    xyz1 = nrlaxis(nrlreverse(tnrls{1}));
    for j=1:2
        rnrlj=nrlreverse(nrls{j});
        xyz = nrlaxis(nrls{j});
        rxyz = nrlaxis(rnrlj);
        if norm(xyz-xyz1)<toler
            tnrls{4}=nrls{j};
            nrls(j) = []; break;
        elseif norm(rxyz-xyz1)<toler
            tnrls{4}=rnrlj;
            nrls(j) = []; break;
        end
    end

    % Find curve 2 whose origin is the end of curve 3
    xyz1 = nrlaxis(nrlreverse(tnrls{3}));
    rnrlj=nrlreverse(nrls{1});
    xyz = nrlaxis(nrls{1});
    if norm(xyz-xyz1)<toler
        tnrls{2}=nrls{1};
    else
        tnrls{2}=rnrlj;
    end
    
end


%% Demo
%  pnts = [ 0.0  3.0  4.5  6.5 8.0 10.0;
%           0.0  0.0  0.0  0.0 0.0  0.0; 
%           2.0  2.0  7.0  4.0 7.0  9.0];   
%  crv1 = nrbmak(pnts, [0 0 0 1/3 0.5 2/3 1 1 1]);
% 
%  pnts= [ 0.0  3.0  5.0  8.0 10.0;
%          10.0 10.0 10.0 10.0 10.0;
%          3.0  5.0  8.0  6.0 10.0];
%  crv2 = nrbmak(pnts, [0 0 0 1/3 2/3 1 1 1]);
% 
%  pnts= [ 0.0 0.0 0.0 0.0;
%          0.0 3.0 8.0 10.0;
%          2.0 0.0 5.0 3.0];
%  crv3 = nrbmak(pnts, [0 0 0 0.5 1 1 1]);
% 
%  pnts= [ 10.0 10.0 10.0 10.0 10.0;
%          0.0   3.0  5.0  8.0 10.0;
%          9.0   7.0  7.0 10.0 10.0];
%  crv4 = nrbmak(pnts, [0 0 0 0.25 0.75 1 1 1]);
% 
%  crv1=nrb2nrl(crv1); crv2=nrb2nrl(crv2); 
%  crv3=nrb2nrl(crv3); crv4=nrb2nrl(crv4);  
%  srf = nrlcoons(crv1, crv2, crv3, crv4);
% 
%  figure; nrlplot(srf,[20 20]);
%  title('Construction of a bilinearly blended Coons surface.');
%  hold off
 
 
 