clear
svgname='hudie2';

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
        R = loadsvg([svgname,'.svg'], 0.1,0 );
        P=R{1,1};
        P(:,2)=max(max(P))-P(:,2);
        Pname=[svgname,'fjy',num2str(length(P)),'.mat'];
end

scatter(P(:,1),P(:,2));

pathname='C:\CodeStore\matlab\vfitcurve\data\';

save ([pathname,Pname],'P')
axis equal