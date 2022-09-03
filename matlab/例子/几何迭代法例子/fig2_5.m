%第二章41页 第三章88页
clear
num=11;%在peaks函数上采样num*num个点
figure
[x,y,z]=peaks(num);
surf(x,y,z)


Pname=['fig2_5-',num2str(num),'.mat'];
pathname='C:\CodeStore\matlab\几何迭代法\data\';

save([pathname,Pname],'x','y','z')

