clear all;
% 启动并行计算
% core_number=2;            %想要调用的处理器个数
% parpool('local',core_number);
tic
example=2;
switch example
    case 1
        % 实验1
        x=0:0.005:1-0.005;
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
    case 3
        %
        x=0:0.005:1-0.005;
        f=1./(0.01+(x-0.3).^2).*(x<0.6)+1./(0.015+(x-0.65).^2).*(x>=0.6);
        n=13;
        Num=length(x);
        a=min(x);
        b=max(x);
        lamda=0.05;
    case 4
        load bird200.mat
        x=gpoint(:,1)';
        f=gpoint(:,2)';
end

f_=f+normrnd(0,1,1,Num);

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
%         ui=[0,0,0,0,0,0,4,4,5,6,6,9.95,9.95,9.95,9.95,9.95,9.95,9.95];


%         u=0:0.005:1-0.005;
%         for j=1:n+1
%             for i=1:length(u)
%                 Njp_u(i,j) = Njp(j, p , u(i), ui);  
%             end
%             subplot(5,5,j),plot(u,Njp_u(:,j));
%         end


        [N,R,P] = kongzhidingdian(M,n,p,x,ui,d);
        plot(x,f_,"*")
        hold on
        DrawSpline(n,p,P,ui,a,b);
        scatter(ui,ones(1,length(ui)))
        [epsilon,e] = shujudianwucha(M,N,P,d);
    case 2
        NP=50;%种群规模
        GM=3;%最大迭代次数
        TournamentSize=3;
        
        maxtime=1;
        for t=1:maxtime
        for i=1:NP
%             X{i}=sort(neijiedian(n,p,a,b)+normrnd(0,0.01,1,n-3));%初始种群并初始化
            pop(i).Position=sort(rand(1,n-3));
%             ui{i}=[zeros(1,p+1) X{i} b*ones(1,p+1)];%方案1 均匀节点向量     
%             [N{i},~,P{i}] = kongzhidingdian(M,n,p,x,ui{i},d);
%             [~,R(i)]=shujudianwucha(M,N{i},P{i},d);
%             BIC(i)=Num*log(1+R(i))+log(Num)*(2*n-p+1);     
            pop(i).Cost=BIC(M,p,x,pop(i).Position,d,b,Num);
        end
        Costs = [pop.Cost];
        [Costs, SortOrder] = sort(Costs);
        pop = pop(SortOrder);
        % Store Best Solution
        BestSol = pop(1);

        % Array to Hold Best Cost Values
        BestCost = zeros(GM, 1);

%         % Store Cost
%         WorstCost = pop(end).Cost;
%         ui_=cell (1,50);
        for gen=1:GM
%             newX=X;

%             轮盘选择操作
            
            fitness_=mean([pop.Cost]);
            minfitness=min([pop.Cost]);
            for m=1:NP/2
                i1 = TournamentSelection(pop, TournamentSize);
                i2 = TournamentSelection(pop, TournamentSize);
                fitness=min([pop(i1).Cost pop(i2).Cost]);
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
                    popc(2*m-1).Position=[pop(i1).Position(pop(i1).Position<cutpoint1),...
                        pop(i2).Position(pop(i2).Position>cutpoint1&pop(i2).Position<cutpoint2),...
                        pop(i1).Position(pop(i1).Position>cutpoint2)];
                    popc(2*m).Position=[pop(i2).Position(pop(i2).Position<cutpoint1),...
                        pop(i1).Position(pop(i1).Position>cutpoint1&pop(i1).Position<cutpoint2),...
                        pop(i2).Position(pop(i2).Position>cutpoint2)];
                else
                    popc(2*m-1).Position=pop(i1).Position;
                    popc(2*m).Position=pop(i2).Position;
                end
%                 n1=length([popc(2*m-1).Position])+p;
%                 n2=length([popc(2*m).Position])+p;
                popc(2*m-1).Cost=BIC(M,p,x,[pop(2*m-1).Position],d,b,Num);
                popc(2*m).Cost=BIC(M,p,x,[pop(2*m).Position],d,b,Num);
            end
%             %交叉操作
%             
%             for m=1:NP/2                
%                 fitness=min([BIC(2*m-1) BIC(2*m)]);                
%                 if fitness<=fitness_
%                     pc=(fitness-minfitness)/(fitness_-minfitness);                    
%                 else
%                     pc=1;
%                 end
%                 if rand()<pc
%                     cutpoint1=a+(b-a)*rand;
%                     cutpoint2=a+(b-a)*rand;
%                     if cutpoint1>cutpoint2
%                         c=cutpoint1;
%                         cutpoint1=cutpoint2;
%                         cutpoint2=c;
%                     end                    
%                     %交叉段
%                     crosstemp1=pop{m*2-1}(pop{m*2-1}>cutpoint1&pop{m*2-1}<cutpoint2);
%                     crosstemp2=pop{m*2}(pop{m*2}>=cutpoint1&pop{m*2}<=cutpoint2);
%                     
%                     temp1=pop{m*2-1}(pop{m*2-1}<cutpoint1|pop{m*2-1}>cutpoint2);
%                     temp2=pop{m*2}(pop{m*2}<=cutpoint1|pop{m*2}>=cutpoint2);
%                     newX{m*2-1}=sort([crosstemp1 temp2]);
%                     newX{m*2}=sort([crosstemp2 temp1]);
% %                     crosstemp1=X{m}(X{m}>cutpoint1&X{m}<cutpoint2);
% %                     crosstemp2=X{NP+1-m}(X{NP+1-m}>=cutpoint1&X{NP+1-m}<=cutpoint2);
% %                     
% %                     temp1=X{m}(X{m}<cutpoint1|X{m}>cutpoint2);
% %                     temp2=X{NP+1-m}(X{NP+1-m}<=cutpoint1|X{NP+1-m}>=cutpoint2);
% %                     newX{m}=sort([crosstemp1 temp2]);
% %                     newX{NP+1-m}=sort([crosstemp2 temp1]);
%                 else 
%                     newX{m*2-1}=pop{m*2-1};
%                     newX{m*2}=pop{m*2};
% %                     newX{m}=X{m};
% %                     newX{NP+1-m}=X{NP+1-m};
%                 end
%             end
            %变异操作
%             for m=1:NP
%                 newui{m}=[zeros(1,p+1) newX{m} b*ones(1,p+1)];%方案1 均匀节点向量  
%                 newn=length(newui{m})-p-2;
%                 [newN{m},~,newP{m}] = kongzhidingdian(M,newn,p,x,newui{m},d);
%                 [~,newR(m)]=shujudianwucha(M,newN{m},newP{m},d);
%                 newBIC(m)=Num*log(1+newR(m))+log(Num)*(2*newn-p+1);                
%             end
            newfitness_=mean([popc.Cost]);
            newminfitness=min([popc.Cost]);
            parfor m=1:NP                
                newfitness=popc(m).Cost;                
                if newfitness<=newfitness_
                    pm=0.5*(newfitness-newminfitness)/(newfitness_-newminfitness);                    
                else
                    pm=0.5;
                end
                iterationcounter=0;
                long=length([popc(m).Position]);
                while iterationcounter<long
                    if rand()<=pm
                        if rand()<=0.5
                            popc(m).Position(end+1)=a+(b-a)*rand;
                            popc(m).Position=sort(popc(m).Position);
                        else
%                             newX{m}(randi(length(newX{m}),1,1))=[];
                            popc(m).Position(randi(length(popc(m).Position),1,1))=[];
                        end
                        iterationcounter=iterationcounter+1;
                    else
                        iterationcounter=iterationcounter+1;
                    end
                end
%                 if length(popc(m).Position)==length(pop(m).Position)
%                     if popc(m).Position==pop(m).Position
%                         popc(m).Position=[];
% %                         popc(m).Cost=[];
% %                     else
% %                         popc(m).Cost=BIC(M,p,x,popc(m).Position,d,b,Num);
%                     end
% %                 else
% %                     popc(m).Cost=BIC(M,p,x,popc(m).Position,d,b,Num);
%                 end
                popc(m).Cost=BIC(M,p,x,popc(m).Position,d,b,Num);
%                 if length(newX{m})==length(pop{m})
%                     if newX{m}==pop{m}
%                         newX{m}=[];
%                     end
%                 end
%                 newn=length(popc(m).Position)+p;
%                 newui{m}=[zeros(1,p+1) newX{m} b*ones(1,p+1)];%方案1 均匀节点向量  
%                 newn=length(newui{m})-p-2;
%                 [newN{m},~,newP{m}] = kongzhidingdian(M,newn,p,x,newui{m},d);
%                 [~,newR(m)]=shujudianwucha(M,newN{m},newP{m},d);
%                 newBIC(m)=Num*log(1+newR(m))+log(Num)*(2*newn-p+1);
                
%                 if newBIC(m)<BIC(m)
%                    X{m}=newX{m};
%                    BIC(m)=newBIC(m);
%                    P{m}=newP{m};
%                    R(m)=newR(m);
%                    ui{m}=newui{m};
%                 end                
            end
            pop=[pop popc];
            Costs = [pop.Cost];
            [Costs, SortOrder] = sort(Costs);
            pop = pop(SortOrder);
            pop = pop(1:NP);
            Costs = Costs(1:NP);
            BestSol = pop(1);
            BestCost(gen) = BestSol.Cost;
            neijiedianshuliang(gen)=length(pop(1).Position);
%             comX=[pop newX];
%             comBIC=[BIC newBIC];
%             comP=[P newP];
%             comR=[R newR];
%             comui=[ui newui];           
%             [sortBIC,index]=sort(comBIC);
%             sortcomX=comX(index);
%             sortcomP=comP(index);
%             sortcomR=comR(index);
%             sortcomui=comui(index);
%             BIC=sortBIC(1:NP);
%             pop=sortcomX(1:NP);
%             P=sortcomP(1:NP);
%             R=sortcomR(1:NP);
%             ui=sortcomui(1:NP);
            
%             [BICbest(gen),position(gen)]=min(BIC);
%             neijiedianshuliang(gen)=length(X{position(gen)});
%             bestP{gen}=P{position(gen)};
%             bestui{gen}=ui{position(gen)};
%             BICbest(gen)=BIC(1);
%             neijiedianshuliang(gen)=length(pop{1});
%             bestP{gen}=P{1};
%             bestui{gen}=ui{1};
        end
            if t==1
                zuishaoneijiedian=neijiedianshuliang(end);
                nbest=1;
                save besttest.mat                
            else
                if neijiedianshuliang(end)<zuishaoneijiedian
                zuishaoneijiedian=neijiedianshuliang(end);
                nbest=1;
                save besttest.mat
                elseif neijiedianshuliang(end)==zuishaoneijiedian
                    nbest=nbest+1;
                    save besttest.mat nbest
                end
            end
            
            
        end
        load besttest.mat
        bestui=[zeros(1,p+1) pop(1).Position b*ones(1,p+1)];
        [~,~,bestP] = kongzhidingdian(M,length(pop(1).Position)+p,p,x,bestui,d);
        toc
        figure
        plot(x,f_,"*")
        hold on
        DrawSpline(neijiedianshuliang(end)+p,p,bestP,bestui,a,b);
        scatter(bestui,ones(1,length(bestui)))
        figure
        plot(BestCost)
        xlabel('迭代次数');
        ylabel('适应度');
        title('适应度进化曲线')
        figure
        plot(neijiedianshuliang)
end
% 关闭并行计算
% delete(gcp('nocreate'));
