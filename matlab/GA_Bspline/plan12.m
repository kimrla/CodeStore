clear all;
example=1;
switch example
    case 1
        % 实验1
        x=0:0.01:1-0.01;
        f=90./(1+exp(-100*(x-0.4)));
        n=8;%控制顶点个数-1
        Num=length(x);%采样点个数
        a=min(x);
        b=max(x);
        lamda=0.025;%节点率
    case 2
        % 实验2
        x=0:0.05:10-0.05;
        f=100./exp(abs(x-5))+(x-5).^2/500;
        n=13;
        Num=length(x);%采样点个数
        a=min(x);
        b=max(x);
        lamda=0.05;%
end

f_=f+normrnd(0,0.01,1,Num);

% plot(x,f_,"*")
% hold on
% plot(x,f)
d=[x;f_]';%给定数据点,di=(xi,yi) i=1~M+1

M=length(x)-1; %论文中是0~M，所以总数length=M+1

p=3;%B样条次数p=3，控制顶点n+1个，节点矢量ui i=1~n+p+2,

plan=2;
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
        plot(x,f_,"*")
        hold on
        DrawSpline(n,p,P,ui,a,b);
        scatter(ui,ones(1,length(ui)))
        [epsilon,e] = shujudianwucha(M,N,P,d);
    case 2
        NP=50;%种群规模
        GM=200;%最大迭代次数
        
        for i=1:NP
%             X{i}=sort(neijiedian(n,p,a,b)+normrnd(0,0.01,1,n-3));%初始种群并初始化
            X{i}=sort(rand(1,n-3));
            ui{i}=[zeros(1,p+1) X{i} b*ones(1,p+1)];%方案1 均匀节点向量     
            [N{i},~,P{i}] = kongzhidingdian(M,n,p,x,ui{i},d);
            [~,R(i)]=shujudianwucha(M,N{i},P{i},d);
            BIC(i)=Num*log(1+R(i))+log(Num)*(2*n-p+1);     
        end
        ui_=cell (1,50);
        for gen=1:GM
%             newX=X;
            %交叉操作
            fitness_=mean(BIC);
            minfitness=min(BIC);
            for m=1:NP/2                
                fitness=min([BIC(2*m-1) BIC(2*m)]);                
                if fitness<=fitness_
                    pc=(fitness-minfitness)/(fitness_-minfitness);                    
                else
                    pc=1;
                end
                if rand()<pc
                    cutpoint1=a+(b-a)*rand;
                    cutpoint2=a+(b-a)*rand;
                    if cutpoint1>cutpoint2
                        c=cutpoint1;
                        cutpoint1=cutpoint2;
                        cutpoint2=c;
                    end                    
                    %交叉段
                    crosstemp1=X{m*2-1}(X{m*2-1}>cutpoint1&X{m*2-1}<cutpoint2);
                    crosstemp2=X{m*2}(X{m*2}>=cutpoint1&X{m*2}<=cutpoint2);
                    
                    temp1=X{m*2-1}(X{m*2-1}<cutpoint1|X{m*2-1}>cutpoint2);
                    temp2=X{m*2}(X{m*2}<=cutpoint1|X{m*2}>=cutpoint2);
                    newX{m*2-1}=sort([crosstemp1 temp2]);
                    newX{m*2}=sort([crosstemp2 temp1]);
                else 
                    newX{m*2-1}=[];
                    newX{m*2}=[];
                end
            end
            %变异操作
            for m=1:NP
                newui{m}=[zeros(1,p+1) newX{m} b*ones(1,p+1)];%方案1 均匀节点向量  
                newn=length(newui{m})-p-2;
                [newN{m},~,newP{m}] = kongzhidingdian(M,newn,p,x,newui{m},d);
                [~,newR(m)]=shujudianwucha(M,newN{m},newP{m},d);
                newBIC(m)=Num*log(1+newR(m))+log(Num)*(2*newn-p+1);                
            end
            newfitness_=mean(newBIC);
            newminfitness=min(newBIC);
            for m=1:NP                
                newfitness=newBIC(m);                
                if newfitness<=newfitness_
                    pm=0.5*(newfitness-newminfitness)/(newfitness_-newminfitness);                    
                else
                    pm=0.5;
                end
                iterationcounter=0;
                long=length(newX{m});
                while iterationcounter<long
                    if rand()<=pm
                        if rand()<=0.5
                            newX{m}(end+1)=a+(b-a)*rand;
                            newX{m}=sort(newX{m});
                        else
                            newX{m}(randi(length(newX{m}),1,1))=[];
                        end
                        iterationcounter=iterationcounter+1;
                    else
                        iterationcounter=iterationcounter+1;
                    end
                end
                newui{m}=[zeros(1,p+1) newX{m} b*ones(1,p+1)];%方案1 均匀节点向量  
                newn=length(newui{m})-p-2;
                [newN{m},~,newP{m}] = kongzhidingdian(M,newn,p,x,newui{m},d);
                [~,newR(m)]=shujudianwucha(M,newN{m},newP{m},d);
                newBIC(m)=Num*log(1+newR(m))+log(Num)*(2*newn-p+1);
%                 if newBIC(m)<BIC(m)
%                    X{m}=newX{m};
%                    BIC(m)=newBIC(m);
%                    P{m}=newP{m};
%                    R(m)=newR(m);
%                    ui{m}=newui{m};
%                 end                
            end
            comX=[X newX];
            comBIC=[BIC newBIC];
            comP=[P newP];
            comR=[R newR];
            comui=[ui newui];           
            [sortBIC,index]=sort(comBIC);
            sortcomX=comX(index);
            sortcomP=comP(index);
            sortcomR=comR(index);
            sortcomui=comui(index);
            BIC=sortBIC(1:NP);
            X=sortcomX(1:NP);
            P=sortcomP(1:NP);
            R=sortcomR(1:NP);
            ui=sortcomui(1:NP);
            
            [BICbest(gen),position(gen)]=min(BIC);
            neijiedianshuliang(gen)=length(X{position(gen)});
            bestP{gen}=P{position(gen)};
            bestui{gen}=ui{position(gen)};
        end
        figure
        plot(x,f_,"*")
        hold on
        DrawSpline(neijiedianshuliang(end)+p,p,bestP{end},bestui{end},a,b);
        scatter(bestui{end},ones(1,length(bestui{end})))
        figure
        plot(BICbest)
        xlabel('迭代次数');
        ylabel('适应度');
        title('适应度进化曲线')
        figure
        plot(neijiedianshuliang)
end
