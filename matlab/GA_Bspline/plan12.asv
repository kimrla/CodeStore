clear all;
example=1;
switch example
    case 1
        % 实验1
        x=0:0.01:1-0.01;
        f=90./(1+exp(-100*(x-0.4)));
        n=8;%控制顶点个数-1
        N=length(x);%采样点个数
        a=min(x);
        b=max(x);
        lamda=0.025;%节点率
    case 2
        % 实验2
        x=0:0.05:10-0.05;
        f=100./exp(abs(x-5))+(x-5).^2/500;
        n=13;
        N=length(x);%采样点个数
        a=min(x);
        b=max(x);
        lamda=0.05;%
end

f_=f+normrnd(0,0.01,1,N);

plot(x,f_,"*")
hold on
% plot(x,f)
d=[x;f_]';%给定数据点,di=(xi,yi) i=1~M+1

M=length(x)-1; %论文中是0~M，所以总数length=M+1

p=3;%B样条次数p=3，控制顶点n+1个，节点矢量ui i=1~n+p+2,

plan=1;
switch plan
    case 1
        ui=jiedianxiangliang(n,p,a,b);%方案1 均匀节点向量


        % u=0:0.01:1;
        % for j=1:n+1
        %     for i=1:length(u)
        %         Njp_u(i,j) = Njp(j, p , u(i), ui);  
        %     end
        %     subplot(5,5,j),plot(u,Njp_u(:,j));
        % end


        [N,R,P] = kongzhidingdian(M,n,p,x,ui,d);

        DrawSpline(n,p,P,ui,a,b);
        scatter(ui,ones(1,length(ui)))
        [epsilon,e] = shujudianwucha(M,N,P,d);
    case 2
        NP=50;%种群规模
        GM=200;%最大迭代次数
        dp=0.6;
        for i=1:NP
            X{i}=sort(neijiedian(n,p,a,b)+normrnd(0,0.01,1,n-3));%初始种群并初始化

            ui{i}=[zeros(1,p+1) X{i} b*ones(1,p+1)];%方案1 均匀节点向量     
            [N{i},~,P{i}] = kongzhidingdian(M,n,p,x,ui{i},d);
            [~,R(i)]=shujudianwucha(M,N{i},P{i},d);
            BIC(i)=Num*log(1+R(i))+log(Num)*(2*n-p+1);     
        end
        ui_=cell (1,50);
        for gen=1:GM
            for m=1:NP
                
            end
        end
end