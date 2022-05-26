load Sgul_wg_kl_ns1.mat thetab n bc 


resfilename=["bimodality_figures_sig48"];
c_theta=thetab;

%set grid points
c_sigmu4=[0:0.05:7];
c_sigmu8=c_sigmu4;

e1=length(c_sigmu4);

sig0val=[0.04,0.5,1.24];

fmat=zeros(e1,e1,length(sig0val));
fmat2=zeros(e1,e1,length(sig0val));


for jj=1:length(sig0val)

c_theta=thetab;
c_theta(26)=sig0val(jj);%loop over cases
    
for j=1:e1
    c_theta(27)=c_sigmu4(j);
    for l=1:e1
    c_theta(28)=c_sigmu8(l);    
    fmat(j,l,jj)=klsgul(thetab,c_theta,n,bc);
    fmat2(j,l,jj)=pfhsgul2(thetab,c_theta,0.05,200,n,bc);
    end
end   

end
save(resfilename)

t1='a)';
t2='b)';
t3='c)';

figure
colormap('parula')

for kk=1:3

    subplot_tight(2,3,kk,0.08);
%subplot_tight(2,2,kk,0.04); %subplot_tight - minimize blank space while maximize figure size
%last input is the desired margin - increasing/decreasing it will increase/decrease blank space
%between subplots and at the edges of the figure
%view(-37.5,60);%-37.5, 30

%set kappa - X, tau - Y:
h=surf(c_sigmu4,c_sigmu8,-fmat(:,:,kk),'EdgeColor','[0.5 0.5 0.5]','EdgeAlpha','0.2'); 

set(h,'LineWidth',0.1)

 %set colormap (default=parula)
%brighten(0.5) %lighten colours (from 0 original, to 1 - 100% lighter)

if kk==1
%title('a $\beta$ =0.90,'Interpreter','latex','FontName','Times New Roman');
%set(gca,'FontName','Times New Roman')
title('a) $\sigma^0_\mu=0.04$','Interpreter','latex');
end


if kk==2
title('b) $\sigma^0_\mu=0.50$','Interpreter','latex');
end

if kk==3
title('c) $\sigma^0_\mu=1.24$','Interpreter','latex');
end

%title([eval(['t' num2str(kk)]) ' \beta = ',num2str(betaval(kk))])

%Xlabel - kappa, ylabel - tau:
ylabel('$\sigma^4_\mu$','FontSize',10,'FontWeight','bold','Interpreter','latex')
xlabel('$\sigma^8_\mu$','FontSize',10,'FontWeight','bold','Interpreter','latex')
xlabh = get(gca,'XLabel');
%set(xlabh,'Position',[12.613 0.956 -3.333])
%set(xlabh,'Position',[12.613 1.556 -3.233])
set(gca,'fontsize',9)


%set(xlabh,'Position',[12.613 0.956 -3.633])
%ylabh = get(gca,'YLabel');
%set(ylabh,'Position',[6.766 1.378 -3.593])
ax=gca;

%caxis(ax, [-ulim dlim]) %set same limits to colorbar
caxis(ax, [-0.015 0])
end

%ax.Position=[0.247 0.029 0.44 0.44]; %adjust position of last subplot to be in the middle
%c=colorbar; %insert colorbar
%c.Position=[0.808 0.029 0.014 0.44]; % adjust position of colorbar
c=colorbar;
c.Position=[0.962 0.532 0.014 0.383]; % adjust position of colorbar a bit in
h=gcf;
set(h,'PaperPositionMode','auto');         
set(h,'PaperOrientation','landscape');
%set(h,'Position',[50 50 1200 800]);
set(h,'Position',[50 50 1000 790]);

set(gcf, 'Color', 'w');
print(gcf, '-dpdf','-r300', 'fig_sig8_kl.pdf','-opengl')



for kk=1:3

    subplot_tight(2,3,kk,0.08);
%subplot_tight(2,2,kk,0.04); %subplot_tight - minimize blank space while maximize figure size
%last input is the desired margin - increasing/decreasing it will increase/decrease blank space
%between subplots and at the edges of the figure
%view(-37.5,60);%-37.5, 30

%set kappa - X, tau - Y:
h=surf(c_sigmu4,c_sigmu8,-fmat2(:,:,kk),'EdgeColor','[0.5 0.5 0.5]','EdgeAlpha','0.2'); 

set(h,'LineWidth',0.1)

 %set colormap (default=parula)
%brighten(0.5) %lighten colours (from 0 original, to 1 - 100% lighter)

if kk==1
%title('a $\beta$ =0.90,'Interpreter','latex','FontName','Times New Roman');
%set(gca,'FontName','Times New Roman')
title('a) $\sigma^0_\mu=0.04$','Interpreter','latex');
end


if kk==2
title('b) $\sigma^0_\mu=0.50$','Interpreter','latex');
end

if kk==3
title('c) $\sigma^0_\mu=1.24$','Interpreter','latex');
end

%title([eval(['t' num2str(kk)]) ' \beta = ',num2str(betaval(kk))])

%Xlabel - kappa, ylabel - tau:
ylabel('$\sigma^4_\mu$','FontSize',10,'FontWeight','bold','Interpreter','latex')
xlabel('$\sigma^8_\mu$','FontSize',10,'FontWeight','bold','Interpreter','latex')
xlabh = get(gca,'XLabel');
%set(xlabh,'Position',[12.613 0.956 -3.333])
%set(xlabh,'Position',[12.613 1.556 -3.233])
set(gca,'fontsize',9)


%set(xlabh,'Position',[12.613 0.956 -3.633])
%ylabh = get(gca,'YLabel');
%set(ylabh,'Position',[6.766 1.378 -3.593])
ax=gca;

%caxis(ax, [-ulim dlim]) %set same limits to colorbar
caxis(ax, [-1 -0.05])
end

%ax.Position=[0.247 0.029 0.44 0.44]; %adjust position of last subplot to be in the middle
%c=colorbar; %insert colorbar
%c.Position=[0.808 0.029 0.014 0.44]; % adjust position of colorbar
c=colorbar;
c.Position=[0.962 0.532 0.014 0.383]; % adjust position of colorbar a bit in
h=gcf;
set(h,'PaperPositionMode','auto');         
set(h,'PaperOrientation','landscape');
%set(h,'Position',[50 50 1200 800]);
set(h,'Position',[50 50 1000 790]);

set(gcf, 'Color', 'w');
print(gcf, '-dpdf','-r300', 'fig_sig8_ed.pdf','-opengl')

save(resfilename)