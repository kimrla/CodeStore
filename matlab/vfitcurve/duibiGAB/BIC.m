function cost = BIC(M,p,x,X,d,b,Num)
%BIC �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
            ui=[zeros(1,p+1) X b*ones(1,p+1)];
            n=length(ui)-p-2;
            [N,~,P] = kongzhidingdian(M,n,p,x,ui,d);
%             [~,R]=shujudianwucha(M,N,P,d);
            R=sum(vecnorm((N*P-d),2,2));
            cost=Num*log(1+R)+log(Num)*(2*n-p+1);
end

