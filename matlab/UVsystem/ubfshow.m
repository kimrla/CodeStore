clear all
x=0:0.001:1;
% for n=1:17
%     for k=1:2^(n-1) 
%         for i=1:length(x)
%             unk_x(i)=unk(n,k,x(i));
%         end
%             subplot,plot(x,unk_x);
%     end
% end
%     

nmax=3;
num=1;
for n=0:nmax
%     figure
    if n==0
     for l=1:length(x)
            unk_x(l)=uni(n,0,x(l));
    end
    subplot(2^nmax,1,1),plot(x,unk_x),ylabel(["n=",n]);
    else
    for i=1:2^(n-1)
        for l=1:length(x)
            unk_x(l)=uni(n,i,x(l));
        end
        num=num+1;
        subplot(2^nmax,1,num),plot(x,unk_x),ylabel(["n=",n,"i=",i]);
    end
    end
end

