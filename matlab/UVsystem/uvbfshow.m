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
n=3;


for k=1:2^(n-1)
    for i=1:length(x)
        unk_x(i)=unk(n,k,x(i));
    end
    subplot(5,5,k),plot(x,unk_x)
end
