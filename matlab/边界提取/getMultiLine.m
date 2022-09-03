clear%%对多条边界取点，给每条边按边界长度占总长比例分配2次幂个点数
svgname='白城市';

qudianjiange=1; %取点模式 0均匀取点 1非均匀取点
switch qudianjiange
    case 0
        pointnum=250; %总取点数目
        R = loadsvg([svgname,'.svg'], 0.001,0 );
        
        P=R{1,1};
        P=[P;P(1,:)];
        % P=unique(R{1,1},'row','stable');
        csl = [0;cumsum(vecnorm(diff(P),2,2))]; %累加长度
        L=csl(end);
        
        t=csl/L;
        tt=linspace(0,1,pointnum)';
        
        Pp=interp1(t,P,tt,'linear');
        clear P
        P=Pp;
        clear Pp
        
        %
        %
        % clear R P
        % R=getpoint(svgname,L,pointnum,0);
        % P=R{1,1};
        
        P(:,2)=max(max(P))-P(:,2);
        Pname=[svgname,'jy',num2str(pointnum),'.mat'];
    case 1
        R = loadsvg([svgname,'.svg'], 0.01,0 );
        %         P=R{1,1};
        for k = 1:length(R)
            Rtemp=unique(R{k},'rows','stable');
            boundary{k} = [Rtemp;R{k}(1,:)];
            
            csl = [0;cumsum(vecnorm(diff(boundary{k}),2,2))];
            L(k)=csl(end);
            t{k}=csl/L(k);
            
            Max(k)=max(max(R{k}));
        end
        High=max(Max);
        LL=cumsum(L);
        nprate=LL/LL(end);%每段边界长度占总长比例的累加值 用于计算每段取点个数
        
        leastn=ceil(log2(length(R)));%2^leastn>=边界个数
        n=7;%预设取点区间为2^n个
        N=max([leastn n]);
        for i=1:length(nprate)
            newnprate(i)=dec2xiaoshu(xiaoshu2dec(nprate(i),N));
        end
        while length(newnprate)-length(unique(newnprate))
            N=N+1;
            for i=1:length(nprate)
                newnprate(i)=dec2xiaoshu(xiaoshu2dec(nprate(i),N));
            end
            
        end
        pointnum=diff([0 newnprate*2^N]);
        
        num=2^N*8;%总点数
        figure
        hold on
        for k = 1:length(R)
            pnum=pointnum(k)*num/2^N;
            %    tt=linspace(0,1,pnum)';
                tt=linspace(1,2*pnum-1,pnum)/(pnum*2)';
                P{k}=interp1(t{k},boundary{k},tt,'linear');
            %     jianju{k}=vecnorm(diff(P{k}),2,2);
            %     canshucha{k}=diff(tt);
%             tt=linspace(0,1,100*length(boundary{k}))';
%             PP=interp1(t{k},boundary{k},tt,'makima');
%             cslPP=[0;cumsum(vecnorm(diff(PP),2,2))];
%             huchang=linspace(1,2*pnum-1,pnum)/(pnum*2)'*cslPP(end);
%             for i=1:length(huchang)
%                 [pMin,minIdx(i)]=min(abs(cslPP-huchang(i)));
%             end
%             P{k}=PP(minIdx,:);
%             
%             hcjianju{k}=diff(huchang);
%             
%             jianju{k}=vecnorm(diff(P{k}),2,2);
%             jianju2{k}=diff(cslPP(minIdx));
%             clear minIdx
            
            P{k}(:,2)=High-P{k}(:,2);
            scatter(P{k}(:,1),P{k}(:,2));
        end
        axis equal
end



pathname='C:\CodeStore\matlab\data\';
Pname=[svgname,'-',num2str(num),'.mat'];
save ([pathname,Pname],'P')

function y=xiaoshu2dec(x,N)
y=zeros(1,N+1);
tempx=x;
for i=1:N+1
    
    y(i)=floor(tempx*2);
    tempx=tempx*2-y(i);
end
if y(end)==1
    y(N)=y(N)+1;
    y(end)=[];
end
for i=fliplr(2:N)
    if y(i)==2
        y(i)=0;
        y(i-1)=y(i-1)+1;
    end
end
% if y(1)==2
%     for i=1:N
%         y(i)=floor(x*2);
%         x=x*2-y(i);
%     end
% end
end

function y=dec2xiaoshu(x)
for i=1:length(x)
    temp(i)=x(i)*2^(-i);
end
y=sum(temp);
end