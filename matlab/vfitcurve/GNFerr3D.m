function GNFerr3D(tzloc,errV,errf,errdb4,errcoif3,errsym5) %分别传入特征点序号，后边的是每种方法每个点的误差
figure
x=linspace(1,length(errV),length(errV));
y=ones(size(errV));
hold on
box on
plot3(x,1*y,errV,'LineWidth',1.1);
plot3(x,2*y,errf,'LineWidth',1.1);
plot3(x,3*y,errdb4,'LineWidth',1.1);
plot3(x,4*y,errcoif3,'LineWidth',1.1);
plot3(x,5*y,errsym5,'LineWidth',1.1);

for i=1:length(tzloc)
    plot3(x(tzloc(i))*y,linspace(0,5,length(errV)),0*y,'k','linewidth',1.5);
end

view(11,17)
xlim([1 length(errV)]);
ylim([0.8 5.2]);
% zlabel('global error')
zlabel('全局误差')
% set(gca,'ytick',[1: 4],'yticklabel',{'V-system','Genetic algorithm','Fourier','Db4 wavelet'},'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
% set(gca,'ytick',[1: 3],'yticklabel',{'V-system','Fourier','Db4 wavelet'},'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
set(gca,'ytick',[1:5],'yticklabel',{'本文方法','Fourier','Db4','Coif3','Sym5'},'linewidth', 1.1, 'fontsize', 10, 'fontname', '微软雅黑')
end