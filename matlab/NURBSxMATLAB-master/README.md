# NURBSxMATLAB

Curve
Input:
u/degree/knotVector/weightVector/controlPoints
Output:
point on curve when @u

Surface
Input:
u/degreeU/knotVectorU
v/degreeV/knotVectorV
weightVector/controlPoints
Output:
point on curve when @u


Curve Input
u: 0~1
degree
knotVector = [0xdegree, ~, 1xdegree];
curve = zeros(length(U),3);
dersCurve = zeros(length(U),degree+1,3);
weightVector = ones(1,length(knotVector)-degree-1);
cvPtX = [0,0,0,0]; % coordinate x
cvPtY = [0,0,0,0];
cvPtZ = [0,0,0,0];
for i = 1:length(U)
u = U(i);
[curve(i,:), dersCurve(i,:,:)] = NurbsCrvPt(u, degree, knotVector, cvPtX, cvPtY, cvPtZ, weightVector);
end
%%
curve(points@u,xyz)
dersCurve(points@u,:,xyz) => derivative from 0 -> degree
curve(points@u,xyz) == dersCurve(points@u,1,xyz)

Surface Input
u: 0~1 (not restrict to 1)
degreeU
knotVectorU = [0xdegree, ~, 1xdegree];
V: 0~1
degreeV
knotVectorV = [0xdegree, ~, 1xdegree];
weightVector = ones(length(knotVectorU)-degreeU-1, length(knotVectorV)-degreeV-1);
surface = zeros(length(U),length(V),3); % 3D pt
dersSurface = zeros(length(U),length(V),degreeU+degreeV+1,degreeU+1,3);
cvPtX = [0,0;0,0]; % coordinate x
cvPtY = [0,0;0,0];
cvPtZ = [0,0;0,0];
for i = 1:length(U)
    u = U(i);
    for j = 1:length(V)
        v = V(j);
        [surface(i,j,:), dersSurface(i,j,:,:,:)] = NurbsSrfPt(u, degreeU, knotVectorU, v, degreeV, knotVectorV, cvPtX, cvPtY, cvPtZ, weightVector);
    end
end
surface(points @u,@v,xyz)
dersSurface(points @u,@v,:,:,xyz) => partial from 0,0 -> degreeU,degreeV
surface(points @u,@v,xyz) == dersSurface(points @u,@v,1,1,xyz)
dersPt( degreeU+1 + degreeV+1 -1, degreeU+1)
