clear
fid=fopen('C:\Users\J\桌面\coord_seligFmt\coord_seligFmt\ag35.dat','r');
dataCell = textscan(fid,'%f %f','headerlines',1);

P = cell2mat(dataCell);
% P=P+1;
P=100*P;
t=canshuhua(P);
k=3;
N=3;
% delta=0.0001;
CList=zeros(2^(N-1)-1,3);
CList(:,1)=1/(2^(N-1)):1/(2^(N-1)):(2^(N-1)-1)/(2^(N-1));
CList(:,2)=CList(:,1);
CList(:,3)=2;
Lambda = LSCurFit_V(P,k,N,t,CList);
Lambdatrapz=LSCurFit_Vtrapz(P,k,N,t,CList);
% Lambda=LSMatrix_V(k,N,t)\P;
Pp=LSMatrix_V(k,N,t)*Lambda;
% wucha=vecnorm((P-Pp),2,2);
% Lambdaupdate=zeros((k+1)*2^(N-2),1);
N=N+1;
NumSeg = 2^(N-2);
for j=1:NumSeg
    [CList1,Lambdaupdate(j:NumSeg:size(Lambda),:)]=VPIA(j,k,N,t,P-Pp,CList);    
end
Lambda=[Lambda;Lambdaupdate];
Aupdate=LSMatrix_V(k,N,t);
Aupdate=Aupdate(:,size(Aupdate,2)/2+1:end);

PpP=Aupdate*Lambdaupdate;

figure,
plot(P(:,1),P(:,2),'.','Color',[255 102 102]/255,'MarkerSize',15);hold on
VCompose(Lambda,k,N)
VCompose(Lambdatrapz,k,N)
function [CList,Lambdaupdate]=VPIA(j,k,N,t,P,CList)
    
    lt=(j-1)/2^(N-2);
    rt=j/2^(N-2);
    
    CList(end+1,:)=[(lt+rt)/2,(lt+rt)/2,2];
    utidx=find((lt<=t)&(t<=rt));
    tupdate=t(utidx);
    Pupdate=P(utidx,:);
    Lambdaupdate=LSCurFit_VPIA(Pupdate,k,N,tupdate,CList,j);
%     Lambdaupdate=LSMatrix_V(k,N,tupdate)\Pupdate;
    CList=CList(end,:);
end
