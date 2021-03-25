
x=0:0.05:10-0.05;
f=100./exp(abs(x-5))+(x-5).^2/500;
g=intergral(f,x-1,x+1,'arrayvalued',true);