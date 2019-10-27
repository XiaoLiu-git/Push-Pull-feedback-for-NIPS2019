function [m] = m_overlapTi(x,xi_c,p_c,Plot)
%m_overlap calculates the overlap between network activities and memory
%pattern xi
%   Plot=0,1
[N,P]=size(xi_c);
x=sign(x*2-1);
xi_c=sign(xi_c*2-1);
m=1/N*(x'*xi_c);% Generat grid 
m=reshape(m,[p_c,P/p_c]);
if Plot==1
figure;
imagesc(m);caxis([-0.1 1 ]);colorbar;set(gca,'linewidth',2);
M=size(m, 1)+1;N=size(m, 2)+1;
hold on;
[xt, yt] = meshgrid(round(linspace(1,M,M)), ...
round(linspace(1, N, N)));%
mesh(yt-0.5, xt-0.5, zeros(size(xt)), 'FaceColor', ...
'None', 'LineWidth', 1, ...
'EdgeColor', 'w');%plot grid
figure_FontSize=20;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',10),'FontSize',figure_FontSize);
title(' M_{1} ')
%saveas(gcf,['/home/lxiao/Downloads/1010/',num2str(ii),'/','moverlap',num2str(ti),'.eps'],'psc2');
end
end