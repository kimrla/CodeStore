%ecg去噪后标记QRS波、R点。
%if 1
%qz=sgolayfilt(M(:,1),2,15);   %ecg去噪 qz去噪后的数据 ，M(1,:)是原始数据
%X=qz; 
%X=m2;
function PEAK_1=getecgtezheng(XC1)
X_1=XC1;  %XC1为信号1去噪后的数据
nOPnum=length(XC1);
for n=2:(nOPnum-1)
    Q0_1(n)=abs(X_1(n+1)-X_1(n-1));    
end
for i=3:(nOPnum-2)
    Q1_1(i)=(Q0_1(i+1)+2*Q0_1(i)+Q0_1(i-1))/4;
end
for j=4:(nOPnum-3)
    Q2_1(j)=abs(X_1(j+2)-2*X_1(j)+X_1(j-2));
    Q3_1(j)=Q1_1(j)+Q2_1(j);
end
ptfs2_1=0.1*max(Q3_1);
stfs2_1=0.0001*max(Q3_1);
FS2_1=[];  %将满足条件的数的下标存入数组FS2中
NFS2_1=[]; %将不满足条件的数的下标存入数组NFS2中
for w=4:(nOPnum-9)
    if (Q3_1(w)>=ptfs2_1&&Q3_1(w+1)>=stfs2_1&&Q3_1(w+2)>=stfs2_1&&Q3_1(w+3)>=stfs2_1&&Q3_1(w+4)>=stfs2_1&&Q3_1(w+5)>=stfs2_1&&Q3_1(w+6)>=stfs2_1)
         FS2_1=[FS2_1,w];
    else
         NFS2_1=[NFS2_1,w];
    end
end
%xhi=find(X>0.2*mean(X));      %采样点
%取顶点_yyp
ataf3_1=0.02;
for n=2:(nOPnum-1)
    M_1(n)=X_1(n+1)-X_1(n-1);
end
AF3_1=[];  %将满足条件的数的下标存入数组AF3中
NAF3_1=[]; %将不满足条件的数的下标存入数组NAF3中
for i=2:(nOPnum-1)
    if (M_1(i)>=ataf3_1&&M_1(i)*X_1(i)>0)
       AF3_1=[AF3_1,i];
    elseif (M_1(i)<=-ataf3_1&&M_1(i)*X_1(i)<0)
       AF3_1=[AF3_1,i];
    else
       NAF3_1=[NAF3_1,i];
    end
end
START_1=[];
for i=1:length(AF3_1)-1
    if((AF3_1(i+1)-AF3_1(i))>2) START_1=[START_1,i];
    else
    end
end
START_1=[START_1 1];
START_1=sort(START_1);
PEAK_1=[];
for i=1:length(START_1)-1
    TEMP_1=[];
    TEMP_1=AF3_1(START_1(i):(START_1(i+1)-1));
    [a,b]=max(X_1(TEMP_1));
    PEAK_1=[PEAK_1,TEMP_1(b)];
end
end