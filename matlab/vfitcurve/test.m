clear
x=0:0.05:10-0.05;
f=100./exp(abs(x-5))+(x-5).^2/500;
y=@(x) 100./exp(abs(x-5))+(x-5).^2/500;
wide=0.6;
for i=1:length(x)
    g(i)=1/wide*integral(y,x(i)-wide/2,x(i)+wide/2);
end
plot(x,f,x,g)



