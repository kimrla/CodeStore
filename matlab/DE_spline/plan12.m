clear all;

tic
example=1;
switch example
    case 1
        % 实验1
        x=0:0.005:1-0.005;
        f=90./(1+exp(-100*(x-0.4)));
        n=8;%控制顶点个数-1
        Num=length(x);%采样点个数
        a=min(x);
        b=max(x);
    case 2
        % 实验2
        x=0:0.05:10-0.05;
        f=100./exp(abs(x-5))+(x-5).^2/500;
        n=13;
        Num=length(x);%采样点个数
        a=min(x);
        b=max(x);
end

f_=f+normrnd(0,1,1,Num);


% 
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


        [Num,R,P] = kongzhidingdian(M,n,p,x,ui,d);

        DrawSpline(n,p,P,ui,a,b);
        scatter(ui,ones(1,length(ui)))
        [epsilon,e] = shujudianwucha(M,Num,P,d);
    case 2
        NP=50;%种群规模NP
        GM=200;%最大迭代次数
        dp=0.6;%删除概率
        maxtime=10;
        
        for t=1:maxtime
            for i=1:NP
                pop(i).position=sort(rand(1,n-3));
                pop(i).cost=BIC(M,p,x,pop(i).position,d,b,Num);
            end
            costs = [pop.cost];
            [costs, sortorder] = sort(costs);
            pop = pop(sortorder);
            % Store Best Solution
            BestSol = pop(1);

            % Array to Hold Best Cost Values
            BestCost = zeros(GM, 1);

    %         % Store Cost
    %         WorstCost = pop(end).Cost;
    %         ui_=cell (1,50);
            for gen=1:GM
                for m=1:NP
                    r=randsample(NP,3);
                    while find(r==m)
                        r=randsample(NP,3);
                    end
                    for i=1:3
                        rand1=rand();
                        npop(i)=pop(r(i));
                        if rand1<dp
                            npop(i).position(randi(length(npop(i).position),1,1))=[];
                        else
                            npop(i).position(end+1)=a+(b-a)*rand;
                            npop(i).position=sort(npop(i).position);
                        end
                        nr(i)=length(npop(i).position);
                    end
                    D=min(nr);
                    for i=1:3
                        for k=1:D
                            tpop(i).position(k)=sum(npop(i).position(k:k+nr(i)-D))/(nr(i)-D+1);
                        end
                    end
                        FG=0.5*exp(GM/(GM+gen)-1);
                        upop(m).position=tpop(1).position+(tpop(2).position-tpop(3).position)*FG;
                        
                        for j=1:D
                            rand4=rand();
                            rand2=rand();
                            rand3=rand();
                            CR=0.5*(rand4+1);
                            if j<=length(pop(m).position) & rand2>CR
                                vpop(m).position(j)=pop(m).position(j);
                            elseif upop(m).position(j)>a &upop(m).position(j)<b
                                vpop(m).position(j)=upop(m).position(j);
                            else
                                vpop(m).position(j)=a+(b-a)*rand3;
                            end
                        end
                        vpop(m).position(j)=sort(vpop(m).position(j));
                        vpop(m).cost=BIC(M,p,x,vpop(m).position,d,b,Num);
                        if vpop(m).cost<pop(m).cost
                            pop(m)=vpop(m);
                        end
                end
%                     pop=[pop vpop];
                    costs = [pop.cost];
                    [costs, sortorder] = sort(costs);
                    pop = pop(sortorder);
%                     pop=pop(1:NP);
%                     costs=costs(1:NP);
                    % Store Best Solution
                    BestSol = pop(1);
                    bestcost(gen)=BestSol.cost;
                    neijiedianshuliang(gen)=length(pop(1).position);
            end
                if t==1
                    zuishaoneijiedian=neijiedianshuliang(end);
                    save besttest.mat
                else
                    if neijiedianshuliang(end)<zuishaoneijiedian
                    zuishaoneijiedian=neijiedianshuliang(end);
                    save besttest.mat
                    end
                end
        end
            load besttest.mat
            bestui=[zeros(1,p+1) pop(1).position b*ones(1,p+1)];
            [~,~,bestP]=kongzhidingdian(M,length(pop(1).position)+p,p,x,bestui,d);
            toc
            figure
            plot(x,f_,"*")
            hold on
            DrawSpline(neijiedianshuliang(end)+p,p,bestP,bestui,a,b);
            scatter(bestui,ones(1,length(bestui)))
            figure
            plot(bestcost)
            xlabel('迭代次数');
            ylabel('适应度');
            title('适应度进化曲线')
            figure
            plot(neijiedianshuliang)
end