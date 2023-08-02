clc;
clearvars;

[p,minTau,tauPure,Tau_vec1,Tau_vec2,x_conv,y_conv]=Run_Minimal_3Tangle_Sequential;
     
%% Plot the results
close all;
Magentas = {[0.4940 0.1840 0.5560],[147,112,219]/255};
Blues    = {[0,104,139]/255,[25,25,112]/255,[65,105,225]/255,[39,64,139]/255,[0,0,128]/255};
Reds     = {[205,55,0]/255,[205,38,38]/255,[139,0,0]/255};
Greens   = {[127,255,212]/255,[102,205,170]/255,[69,139,116]/255}  ;

fntsize  = 20;
rows     = 2;
cols     = 1;

subplot(rows,cols,1)

plot(p,tauPure,'marker','o','color',Magentas{1},...
    'markerfacecolor',Magentas{1},'markersize',8)
hold on
plot(p,Tau_vec1,'marker','^','color',Blues{3},'markerfacecolor',Blues{3},'markersize',6)
% hold on
% plot(p,Tau_vec2,'marker','d')

ylim([min(Tau_vec1)-0.001,1])
fig_props(fntsize,'p','$\tau_3$ (Pure)')

legend({'$\tau_3(|\psi\rangle)$','$\tau_3(|v_\pm\rangle)$','$\tau_3(|v_-\rangle)$'},...
       'interpreter','latex','location','best','color','none','edgecolor','none')

set(gca,'ytick',[0.98,0.99,1],'yticklabels',{'0.98','0.99','1'})
xlim([min(p),max(p)])
xlabel('$p$','interpreter','latex')
%============================================================================
subplot(rows,cols,2)

plot(p,minTau,'marker','o','markersize',8,'markerfacecolor',Reds{1},'color',Reds{1})
hold on
plot(x_conv,y_conv,'linewidth',1,'linestyle','--','color','k','marker','*')

fig_props(fntsize,'p','$\tau_3$ (Mixed)')
xlabel('$p$','interpreter','latex')
legend({'min$_\chi(\tau_3)$','convex hull'},...
       'interpreter','latex','location','best','color','none','edgecolor','none')
 
ylim([min(minTau)-0.001,1])
xlim([min(p),max(p)])








