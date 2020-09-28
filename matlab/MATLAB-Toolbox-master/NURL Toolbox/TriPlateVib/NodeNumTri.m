function [Num, BN, TN]=NodeNumTri(num1, num2)

TNe=num1*(num2-1)+1;
Num1=reshape(1:num1*num2, [num1, num2]);
Num1(:,end)=Num1(1,end);
Num2=0*Num1; Num3=0*Num1; Num4=0*Num1; 
Num2(:,end)=Num1(1,end);
Num3(:,end)=Num1(1,end);
Num4(:,end)=Num1(1,end);
Num2(1,:)=Num1(end,:);
N1=Num1(1,end);
Num2(2:end, 1:end-1)=reshape(N1+(1:(num1-1)*(num2-1)), [num1-1, num2-1]);
N2=Num2(end, end-1);
Num3(1,:)=Num2(end,:);
Num3(2:end, 1:end-1)=reshape(N2+(1:(num1-1)*(num2-1)), [num1-1, num2-1]);
N3=Num3(end, end-1);
Num4(1,:)=Num3(end,:);
Num4(end,:)=Num1(1,:);
Num4(2:end-1, 1:end-1)=reshape(N3+(1:(num1-2)*(num2-1)), [num1-2, num2-1]);
BN=[Num1(:,1); Num2(:,1); Num3(:,1); Num4(:,1)];
BN=RemDuplicate(BN);
TN=Num4(end-1, end-1);
Num1=Num1(1:TNe);
Num2=Num2(1:TNe);
Num3=Num3(1:TNe);
Num4=Num4(1:TNe);

Num=[Num1; Num2; Num3; Num4];




