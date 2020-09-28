clear
clc

%% varify BasisFuns
% degree = 2;
% knotVector = [0,0,0,1,2,3,4,4,5,5,5];
% u = 5/2;
% [basisFunValue, dersBasisFunValue] = BasisFuns(u, degree, knotVector);

%% Surface
% degreeU = 2;
% knotVectorU = [0, 0, 0, 7.7060181011129467, 15.412036202225893, 15.412036202225893, 15.412036202225893];
% U = linspace(0,15.412036202225893);
% degreeV = 2;
% knotVectorV = [0, 0, 0, 7.1181589352563979, 14.236317870512796, 14.236317870512796, 14.236317870512796];
% V = linspace(0,14.236317870512796);
% weightVector = ones(length(knotVectorU)-degreeU-1, length(knotVectorV)-degreeV-1);
% cvPtX = [ -5.6154732922582369, 11.300198435115878, 14.643115727078438, 18.634377426229133;
%     6.7902280829194961, 9.0433806443932259, 10.730356052276349, 10.025771485038982;
%     2.3684705522552045, 3.9477120476917031, 13.238982721975979, 16.854983730304763;
%     1.8753182909394326, -2.2959743902089347, 8.4852498635469118, 8.7736178845056507];
% cvPtY = [ -2.3406130349086478, 5.883021828144388, -2.2304178664224148, -2.1713296092996632;
%     11.520420269598095, 4.5890772299414895, 8.2998590678519939, 11.021817120658726;
%     15.803463281522482, 18.479549414893857, 16.764695944660776, 18.58155205594046;
%     17.087864131508812, 24.284003784872557, 14.937881529902, 13.021857839106911];
% cvPtZ = [ 0, 0, 0, 0;
%     -0.45358741447932982, -0.90149442591237283, 0.15361649200790986, 0;
%     1.789899339917453, 3.5573832658739941, -0.60618537655198212, 0;
%     3.8351814758997342, 7.6223338931366005, -1.2988612684893679, 0];
% surface = zeros(length(U),length(V),3); % 3D pt
% dersSurface = zeros(length(U),length(V),degreeU+degreeV+1,degreeU+1,3);
% for i = 1:length(U)
%     u = U(i);
%     for j = 1:length(V)
%         v = V(j);
%         [surface(i,j,:), dersSurface(i,j,:,:,:)] = NurbsSrfPt(u, degreeU, knotVectorU, v, degreeV, knotVectorV, cvPtX, cvPtY, cvPtZ, weightVector);
%     end
% end
% xSc = reshape(surface(:,:,1),[],1);
% ySc = reshape(surface(:,:,2),[],1);
% zSc = reshape(surface(:,:,3),[],1);
% scatter3(xSc, ySc, zSc);
% hold on
% scatter3(reshape(cvPtX,[],1), reshape(cvPtY,[],1), reshape(cvPtZ,[],1));
% 
% fileID = fopen('SurfaceLisp.lsp','w');
% for i = 1:length(U)
%     for j = 1:length(V)
%         fprintf(fileID,'(COMMAND "POINT" "%f,%f,%f")\n',surface(i,j,1),surface(i,j,2),surface(i,j,3));
%     end
% end
% fclose('all');

%% Surface test
% degreeU = 2;
% knotVectorU = [0, 0, 0, 1/2, 1, 1, 1];
% U = linspace(0,1);
% degreeV = 2;
% knotVectorV = [0, 0, 0, 1/2, 1, 1, 1];
% V = linspace(0,1);
% weightVector = ones(length(knotVectorU)-degreeU-1, length(knotVectorV)-degreeV-1);
% cvPtX = [ 1, 1, 1, 1;
%     2, 2, 2, 2;
%     3, 3, 3, 3;
%     4, 4, 4, 4];
% cvPtY = [ 1, 2, 3, 4;
%     1, 2, 3, 4;
%     1, 2, 3, 4;
%     1, 2, 3, 4];
% cvPtZ = [ 1, 1, 1, 1;
%     1, 1, 1, 1;
%     2, 2, 2, 2;
%     2, 2, 2, 2];
% surface = zeros(length(U),length(V),3); % 3D pt
% for i = 1:length(U)
%     u = U(i);
%     for j = 1:length(V)
%         v = V(j);
%         surface(i,j,:) = NurbsSrfPt(u, degreeU, knotVectorU, v, degreeV, knotVectorV, cvPtX, cvPtY, cvPtZ, weightVector);
%     end
% end
% xSc = reshape(surface(:,:,1),[],1);
% ySc = reshape(surface(:,:,2),[],1);
% zSc = reshape(surface(:,:,3),[],1);
% scatter3(xSc, ySc, zSc);
% hold on
% scatter3(reshape(cvPtX,[],1), reshape(cvPtY,[],1), reshape(cvPtZ,[],1));

%% Curve
U = linspace(0,1);
degree = 3;
curve = zeros(length(U),3); % 3D pt
dersCurve = zeros(length(U),degree+1,3);
knotVector = [0, 0, 0, 0, 1, 1, 1, 1];
weightVector = ones(1,length(knotVector)-degree-1);
cvPtX = [0.084641712804496352,...
    3.3218617663495729,...
    6.9286011411941901,...
    2.3699679585258182];
cvPtY = [0.063481284603371504,...
    5.0654020277678935,...
    -3.0670589716973753,...
    -3.9570000736100668];
cvPtZ = [0, 0, 0, 0];
cvPt = [cvPtX', cvPtY', cvPtZ'];
for i = 1:length(U)
u = U(i);
[curve(i,:), dersCurve(i,:,:)] = NurbsCrvPt(u, degree, knotVector, cvPt, weightVector);
end
plot3(curve(:,1),curve(:,2),curve(:,3));
hold on
scatter3(cvPtX, cvPtY, cvPtZ);

fileID = fopen('CurveLisp.lsp','w');
for i = 1:length(U)
    fprintf(fileID,'(COMMAND "POINT" "%f,%f")\n',curve(i,1),curve(i,2));
end
fclose('all');
