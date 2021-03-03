example=1;
switch example
    case 1
        % 实验1
        x=0:0.005:1-0.005;
        f=90./(1+exp(-100*(x-0.4)));        
        Num=length(x);%采样点个数
        
    case 2
        % 实验2
        x=0:0.05:10-0.05;
        f=100./exp(abs(x-5))+(x-5).^2/500;        
        Num=length(x);%采样点个数  
        
    case 3
        %
        x=0:0.005:1-0.005;
        f=1./(0.01+(x-0.3).^2).*(x<0.6)+1./(0.015+(x-0.65).^2).*(x>=0.6);
        Num=length(x);
        
end
f_=f+normrnd(0,1,1,Num);
d=[x;f_]';