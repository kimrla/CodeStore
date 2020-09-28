%% ========================================================================
function polygon = IntersectionPolygon(edgeMat, points, dIdx, debug)
% edgeMat is an edge intersection matrix with 3 rows for edges between
% the points 1-3, 1-2, & 2-3 of the triangle 1 and 3 columns for the same
% edges of the triangle 2. If 2 edges intersect a point of intersection
% is calculated and stored in array "points" followed by points of the
% triangles 1 & 2.  This function calculates the polygon of the intersection
% between 2 triangles.

persistent orderLUT verified
if isempty(orderLUT) || isempty(orderLUT{3})
  % This pre-calculated look-up table is used to quickly look up the order of
  % the vertices in array "points" which make up polygon of the intersection
  % between 2 triangles. A unique key is calculated for each edgeMat using
  % dot product between edgeMat(:) and [256 128 64 32 16 8 4 2 1], which is
  % used to look up point order around the polygon. Negative numbers in the
  % LUT indicate values which were not observed yet so they were not
  % independently verified.
  % reshape(sprintf('%09s',dec2base(key, 2)),3,3) will convert from the key
  % to matrix.
  OrderLUT = zeros(432,1);  
  OrderLUT(003) = 127;
  OrderLUT(005) = 128;
  OrderLUT(006) = 126;
  OrderLUT(009) = 124;
  OrderLUT(010) = 1427;
  OrderLUT(012) = 1428;
  OrderLUT(017) = 1427;
  OrderLUT(018) = 124;
  OrderLUT(020) = 1426;
  OrderLUT(024) = 127;
  OrderLUT(027) = 1243;
  OrderLUT(029) = 12438;
  OrderLUT(030) = 12034;
  OrderLUT(033) = 1428;
  OrderLUT(034) = 1426;
  OrderLUT(036) = 124;
  OrderLUT(040) = 128;
  OrderLUT(043) = 21834;
  OrderLUT(045) = 1243;
  OrderLUT(046) = 21349;
  OrderLUT(048) = 126;
  OrderLUT(051) = 12340;
  OrderLUT(053) = 12943;
  OrderLUT(054) = 1243;
  OrderLUT(065) = 125;
  OrderLUT(066) = 1527;
  OrderLUT(068) = 1825;
  OrderLUT(072) = 123;
  OrderLUT(080) = 1327;
  OrderLUT(083) = 15234;
  OrderLUT(085) = -15234;
  OrderLUT(086) = -15243;
  OrderLUT(090) = 13247;
  OrderLUT(092) = -13247;
  OrderLUT(096) = 1328;
  OrderLUT(099) = 152834;
  OrderLUT(101) = 15234;
  OrderLUT(102) = 152349;
  OrderLUT(106) = 132847;
  OrderLUT(108) = 13247;
  OrderLUT(114) = 102347;
  OrderLUT(116) = -13247;
  OrderLUT(129) = 1527;
  OrderLUT(130) = 125;
  OrderLUT(132) = 1526;
  OrderLUT(136) = 1327;
  OrderLUT(139) = 15243;
  OrderLUT(141) = 152438;
  OrderLUT(142) = 152034;
  OrderLUT(144) = 123;
  OrderLUT(153) = 12347;
  OrderLUT(156) = 123047;
  OrderLUT(160) = 1326;
  OrderLUT(163) = -152043;
  OrderLUT(165) = 13247;
  OrderLUT(166) = 15234;
  OrderLUT(169) = -182347;
  OrderLUT(172) = 193247;
  OrderLUT(177) = -132047;
  OrderLUT(180) = 13247;
  OrderLUT(192) = 127;
  OrderLUT(195) = 1243;
  OrderLUT(197) = 12438;
  OrderLUT(198) = 12034;
  OrderLUT(202) = 12364;
  OrderLUT(204) = 123648;
  OrderLUT(209) = 21364;
  OrderLUT(212) = -21364;
  OrderLUT(216) = 1243;
  OrderLUT(225) = -124638;
  OrderLUT(226) = 120364;
  OrderLUT(232) = 12438;
  OrderLUT(238) = 124356;
  OrderLUT(240) = 12034;
  OrderLUT(245) = -214356;
  OrderLUT(257) = 1528;
  OrderLUT(258) = 1526;
  OrderLUT(260) = 125;
  OrderLUT(264) = 1328;
  OrderLUT(267) = -152438;
  OrderLUT(269) = 15243;
  OrderLUT(270) = -152943;
  OrderLUT(272) = 1326;
  OrderLUT(275) = 152340;
  OrderLUT(277) = 152943;
  OrderLUT(278) = 15243;
  OrderLUT(281) = 182347;
  OrderLUT(282) = -103247;
  OrderLUT(288) = 123;
  OrderLUT(297) = 12347;
  OrderLUT(298) = -123947;
  OrderLUT(305) = 123947;
  OrderLUT(306) = 12347;
  OrderLUT(320) = 128;
  OrderLUT(323) = 21834;
  OrderLUT(325) = 1243;
  OrderLUT(326) = 21349;
  OrderLUT(330) = -123648;
  OrderLUT(332) = 12364;
  OrderLUT(337) = 183642;
  OrderLUT(340) = -129364;
  OrderLUT(344) = 21834;
  OrderLUT(350) = -124365;
  OrderLUT(353) = 12463;
  OrderLUT(354) = 136492;
  OrderLUT(360) = 1243;
  OrderLUT(368) = 12943;
  OrderLUT(371) = 126543;
  OrderLUT(384) = 126;
  OrderLUT(387) = 12340;
  OrderLUT(389) = 12943;
  OrderLUT(390) = 1243;
  OrderLUT(394) = -103642;
  OrderLUT(396) = 129364;
  OrderLUT(401) = 123640;
  OrderLUT(404) = 12364;
  OrderLUT(408) = 12340;
  OrderLUT(413) = 215643;
  OrderLUT(417) = -136492;
  OrderLUT(418) = 12463;
  OrderLUT(424) = 13492;
  OrderLUT(427) = -213456;
  OrderLUT(432) = 1342;
  
  % Convert to more convinient format
  orderLUT = cell(size(OrderLUT));
  for i = 1:size(OrderLUT,1)
    polygon = abs(OrderLUT(i));
    if polygon>0
      polygon = num2str(polygon)-48; % Convert from a single number to array of digits
      polygon(polygon==0) = 10;      % 0 stands for 10
      orderLUT{i} = polygon;
    end
  end
  % Negative numbers in the LUT indicate values which were not observed yet
  % so they were not independently verified.
  verified = OrderLUT>0;
  clear OrderLUT
end

%% Calculate unique key for each edgeMat configuration
key = dot(1*edgeMat(:)', [256 128 64 32 16 8 4 2 1]);
assert(key<=432, 'Error: in IntersectionPolygon: key is out of bound');

%% Look up the point order around the polygon
polygon = orderLUT{key};
if (isempty(polygon))
  return
end

%% in a rare case of 2 intersections there is ambiguity if one or two
% vertices of the triangle lay inside the other triangle. OrderLUT stores
% only the single vertex cases.
nx = nnz(edgeMat(:));
if nx==2
  pList = polygon;       % list of vertices to check
  pList(pList<=nx) = []; % keep only the triangle points of the polygon
  flip = false;    % was there a flip from single vertex to vertices case?
  for ip = 1:length(pList)
    p = pList(ip);                 % point to check
    t = floor((p-nx-1)/3);         % does it belong to triangle 0 or 1 (actually 1 or 2)
    tri = (1:3) + nx + 3*abs(1-t); % Points belonging to the other triangle
    if ~PointInTriangle2D(points(p,dIdx), points(tri,dIdx))
      d = nx+t*3;    % offset
      % "p-d" is vertex number of point just tested: 1, 2, or 3. "b, c" are
      % the other 2 vertices
      [b, c] = otherDim(p-d);
      polygon = [polygon(polygon~=p), b+d, c+d]; % remove i2 and add i0 and i1
      flip = true;
    end
  end
  if flip
    % if ther were any flips than use existing codes to figure out the
    % order of the points around the polygon
    DT = delaunayTriangulation(points(polygon,dIdx));
    idx = freeBoundary(DT)';
    idx(2,:) = [];
    polygon = polygon(idx);
  end
end

%% Check to duplicate points
tol = 1e6;
P = round(points(polygon,:)*tol)/tol;
[~,ia] = unique(P,'rows'); % V = P(ia,:) and P = V(ic,:).
polygon = polygon(sort(ia));

%% Test the results using more expensive function
doPlot = (~verified(key));
if debug && length(polygon)>3
  DT = delaunayTriangulation(points(polygon,dIdx));
  idx = freeBoundary(DT)';
  idx(2,:) = [];
  k = max(abs(diff(idx)));
  %doPlot = (k>1 && k<(length(idx)-1)) || (~verified(key));
  assert(k==1 || k==(length(idx)-1), 'Two triangle intersection polygon is not convex')
end
if debug && doPlot % plot the interesting cases
  PlotTwoTriangles(points, polygon, 'm')
  title(sprintf('key = %i', key));
end 

end % function