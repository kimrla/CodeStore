%% Compute point on non-uniform rational B-spline surface
% Input:  u,v, degree(U,V), knotVector(U,V), cvPt(X,Y,Z), weightVector
% Output: point on non-uniform rational B-spline surface
% column:u / row:v
% it doesn't work at boundaries.
% Xu Yi, 2019

%%
function [surfacePt, dersSurfacePt] = NurbsSrfPt(u, degreeU, knotVectorU, v, degreeV, knotVectorV, cvPtX, cvPtY, cvPtZ, weightVector)
% u
knotspanIndexU = FindSpan(u, knotVectorU);
[basisFunValueU, dersBasisFunValueU] = BasisFuns(u, degreeU, knotVectorU);
% v
knotspanIndexV = FindSpan(v, knotVectorV);
[basisFunValueV, dersBasisFunValueV] = BasisFuns(v, degreeV, knotVectorV);
% w.*P
WcvPtX = cvPtX.* weightVector;
WcvPtY = cvPtY.* weightVector;
WcvPtZ = cvPtZ.* weightVector;
if knotspanIndexU == length(knotVectorU) - degreeU && knotspanIndexV == length(knotVectorV) - degreeV % special case
    PtX = cvPtX(end,end);
    PtY = cvPtY(end,end);
    PtZ = cvPtZ(end,end);
elseif knotspanIndexU == length(knotVectorU) - degreeU % special case / like curve @ V
    W = basisFunValueV(end,:)...
        * ( weightVector( end , (knotspanIndexV-degreeV):knotspanIndexV )' );
    PtX = basisFunValueV(end,:)...
        * ( WcvPtX( end , (knotspanIndexV-degreeV):knotspanIndexV )' )...
        / W;
    PtY = basisFunValueV(end,:)...
        * ( WcvPtY( end , (knotspanIndexV-degreeV):knotspanIndexV )' )...
        / W;
    PtZ = basisFunValueV(end,:)...
        * ( WcvPtZ( end , (knotspanIndexV-degreeV):knotspanIndexV )' )...
        / W;
elseif knotspanIndexV == length(knotVectorV) - degreeV % special case / like curve @ U
    W = basisFunValueU(end,:)...
        * ( weightVector( (knotspanIndexU-degreeU):knotspanIndexU , end ) );
    PtX = basisFunValueU(end,:)...
        * ( WcvPtX( (knotspanIndexU-degreeU):knotspanIndexU , end ) )...
        / W;
    PtY = basisFunValueU(end,:)...
        * ( WcvPtY( (knotspanIndexU-degreeU):knotspanIndexU , end ) )...
        / W;
    PtZ = basisFunValueU(end,:)...
        * ( WcvPtZ( (knotspanIndexU-degreeU):knotspanIndexU , end ) )...
        / W;
else
    W = basisFunValueU(end,:)...
        * weightVector( (knotspanIndexU-degreeU):knotspanIndexU, (knotspanIndexV-degreeV):knotspanIndexV )...
        * basisFunValueV(end,:)';
    PtX = basisFunValueU(end,:)...
        * WcvPtX( (knotspanIndexU-degreeU):knotspanIndexU, (knotspanIndexV-degreeV):knotspanIndexV )...
        * basisFunValueV(end,:)'...
        / W;
    PtY = basisFunValueU(end,:)...
        * WcvPtY( (knotspanIndexU-degreeU):knotspanIndexU, (knotspanIndexV-degreeV):knotspanIndexV )...
        * basisFunValueV(end,:)'...
        / W;
    PtZ = basisFunValueU(end,:)...
        * WcvPtZ( (knotspanIndexU-degreeU):knotspanIndexU, (knotspanIndexV-degreeV):knotspanIndexV )...
        * basisFunValueV(end,:)'...
        / W;
end
surfacePt = [PtX, PtY, PtZ];
% derivative. The NURBS book 2nd page136 (4.20)
% initialize / combination
dersBasisCombU = zeros(degreeU+1,degreeU+1);
dersBasisCombU(1,:) = basisFunValueU(end,:);
dersBasisCombU(2:end,:) = dersBasisFunValueU;
dersBasisCombV = zeros(degreeV+1,degreeV+1);
dersBasisCombV(1,:) = basisFunValueV(end,:);
dersBasisCombV(2:end,:) = dersBasisFunValueV;
% boundary
if knotspanIndexU == length(knotVectorU) - degreeU % special case / like curve @ V
    [knotspanIndexU, dersBasisCombU] = specialCase(knotspanIndexU, dersBasisCombU, degreeU);
end
if knotspanIndexV == length(knotVectorV) - degreeV % special case / like curve @ U
    [knotspanIndexV, dersBasisCombV] = specialCase(knotspanIndexV, dersBasisCombV, degreeV);
end
% Aders
AdersX = genAders(WcvPtX,...
    degreeU, knotspanIndexU, dersBasisCombU,...
    degreeV, knotspanIndexV, dersBasisCombV);
AdersY = genAders(WcvPtY,...
    degreeU, knotspanIndexU, dersBasisCombU,...
    degreeV, knotspanIndexV, dersBasisCombV);
AdersZ = genAders(WcvPtZ,...
    degreeU, knotspanIndexU, dersBasisCombU,...
    degreeV, knotspanIndexV, dersBasisCombV);
% wders
wders = genWders(weightVector,...
    degreeU, knotspanIndexU, dersBasisCombU,...
    degreeV, knotspanIndexV, dersBasisCombV);
% dersSurfacePt
dersPtX = genDersPt(AdersX, degreeU, degreeV, wders);
dersPtY = genDersPt(AdersY, degreeU, degreeV, wders);
dersPtZ = genDersPt(AdersZ, degreeU, degreeV, wders);
dersSurfacePt = zeros(degreeU+degreeV+1,degreeU+1,3);
dersSurfacePt(:,:,1) = dersPtX;
dersSurfacePt(:,:,2) = dersPtY;
dersSurfacePt(:,:,3) = dersPtZ;
end

function [knotspanIndex, dersBasisComb] = specialCase(knotspanIndex, dersBasisComb, degree)
knotspanIndex = knotspanIndex -1;
dersBasisCombU_temp = dersBasisComb;
dersBasisComb(:,1) = zeros(degree+1,1);
dersBasisComb(:,2:end) = dersBasisCombU_temp(:,1:end-1);
end
function Aders = genAders(WcvPt,...
    degreeU, knotspanIndexU, dersBasisCombU,...
    degreeV, knotspanIndexV, dersBasisCombV)
Aders = zeros(degreeU+degreeV+1,degreeU+1);
for i = 1:degreeU+1
    for j = 1:degreeV+1
        Aders(i+j-1,i) = dersBasisCombU(i,:)... % index from i.
            * WcvPt( (knotspanIndexU-degreeU):knotspanIndexU, (knotspanIndexV-degreeV):knotspanIndexV )...
            * dersBasisCombV(j,:)';
    end
end
end
function wders = genWders(weightVector,...
    degreeU, knotspanIndexU, dersBasisCombU,...
    degreeV, knotspanIndexV, dersBasisCombV)
wders = zeros(degreeU+degreeV+1,degreeU+1);
for i = 1:degreeU+1
    for j = 1:degreeV+1
        wders(i+j-1,i) = dersBasisCombU(i,:)... % index from i.
            * weightVector( (knotspanIndexU-degreeU):knotspanIndexU, (knotspanIndexV-degreeV):knotspanIndexV )...
            * dersBasisCombV(j,:)';
    end
end
end
function dersPt = genDersPt(Aders, degreeU, degreeV, wders)
dersPt = zeros(degreeU+degreeV+1,degreeU+1);
for i = 1:degreeU+1
    for j = 1:degreeV+1
        wPt = 0; % initialize
        if i == 1 % special case
        elseif j == 1
        else
            for ii = 2:i %
                wPt_temp = factorial(i-1)/factorial(ii-1)/factorial(i-ii)...
                    * wders(ii,ii) * dersPt(i-ii+j,i-ii+1);
                wPt = wPt + wPt_temp;
            end
        end
        if j == 1 % special case
        elseif i == 1
        else
            for jj = 2:j
                wPt_temp = factorial(j-1)/factorial(jj-1)/factorial(j-jj)...
                    * wders(jj,1) * dersPt(i+j-jj,i);
                wPt = wPt + wPt_temp;
            end
        end
        if i == 1 % special case
        elseif j == 1 % special case
        else
            for ii = 2:i
                for jj = 2:j
                    wPt_temp = factorial(i-1)/factorial(ii-1)/factorial(i-ii)...
                        * factorial(j-1)/factorial(jj-1)/factorial(j-jj)...
                        * wders(ii+jj-1,ii) * dersPt(i-ii+1+j-jj+1-1,i-ii+1);
                    wPt = wPt + wPt_temp;
                end
            end
        end
        dersPt(i+j-1,i) = ( Aders(i+j-1,i) - wPt ) / wders(1,1);
    end
end
end
