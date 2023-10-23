%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------

clear
clc

s0=0;
s1=-1;

X=[0 1; 1 0];
Y=[0 -1i ; 1i 0];
Z=[1 0 ; 0 -1];
Ze=[s0 0 ; 0 s1];
Rx=@ (err) expm(-1i*(pi+err)/2*X);
Ry=@ (err) expm(-1i*(pi+err)/2*Y);

%Average the fidelity over arbitrary initial states to find the gate fid of
%the pi-pulse

ket0 = [1;0];
ket1 = [0;1];
ket2 = 1/sqrt(2)*[1;1];
ket3 = 1/sqrt(2)*[1;-1];
ket4 = 1/sqrt(2)*[1;1i];
ket5 = 1/sqrt(2)*[1;-1i];

kets = {ket0,ket1,ket2,ket3,ket4,ket5};

for jj=1:length(kets)
   
    rho0s{jj}=kets{jj}*kets{jj}';
    
end


%% Systematic pulse error

c = 2*pi*1e-3;

B0        = 403;
B0        = B0*1e-4;   
t         = 3; %\mu s
gamma_e   = 1.760859*10^11 ;  %rad/T*s

we = gamma_e * B0 ;               % Larmor freq of electron
we = we*1e-6;                       %  MHz
we = we*1e3;                        % in KHz
we = we/(2*pi);                     

we       = c*we;
iter_Max = 300;

pulse_error1 = 2/100;
pulse_error2 = 8/100;

Hfree = we*Ze;
Uf    = @(t) expm(-1i*t*Hfree); 

UCPMG=@(t,err) Uf(t/4)*Rx(err(1))*Uf(t/2)*Rx(err(2))*Uf(t/4);
UXY2=@(t,err)  Uf(t/4)*Rx(err(1))*Uf(t/2)*Ry(err(2))*Uf(t/4);

for N=1:iter_Max
    
    Uideal_CPMG = UCPMG(t,[0,0])^N;
    Uideal_XY2  = UXY2(t,[0,0])^N;
    
    Uer_CPMG    = UCPMG(t,[pulse_error1,pulse_error1])^N;
    Uer_XY2     = UXY2(t,[pulse_error1,pulse_error1])^N;
    Uer_XY2_V2  = UXY2(t,[pulse_error2,pulse_error2])^N;
    
    for ll=1:length(rho0s)
       
        r0 = rho0s{ll};
        
        rj_CPMG    = Uer_CPMG*r0*Uer_CPMG';
        rj_XY2     = Uer_XY2*r0*Uer_XY2';
        rj_XY2_V2  = Uer_XY2_V2*r0*Uer_XY2_V2';
        
        fCPMG(ll)   = real(trace(Uideal_CPMG*r0*Uideal_CPMG'*rj_CPMG));
        fXY2(ll)    = real(trace(Uideal_XY2*r0*Uideal_XY2'*rj_XY2));
        fXY2_V2(ll) = real(trace(Uideal_XY2*r0*Uideal_XY2'*rj_XY2_V2));
    end
    
    FidCPMG(N)    = sum(fCPMG)/6;
    FidXY2(N)     = sum(fXY2)/6;
    FidXY2_V2(N)  = sum(fXY2_V2)/6;
    
end
close all

h=semilogy(1:iter_Max,[1-FidXY2;1-FidXY2_V2],'linewidth',2)
h(1).Marker='o';
%ylim([0.95,1])

legend('XY2 $\epsilon=2\%$','XY2 $\epsilon=8\%$',...
       'location','best','color','none',...
       'edgecolor','none','interpreter','latex')
xlabel('$N$','interpreter','latex')
ylabel('$1-F$','interpreter','latex')

set(gca,'fontsize',24,'fontname','Microsoft Sans Serif')

set(gcf,'color','w')

set(gca, 'Box', 'on', 'TickDir', 'out', 'TickLength', [.003 .003], ...
    'XMinorTick', 'off', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
    'LineWidth', 0.5)



%% Random pulse errors
clear
clc
reps = 1000;

Nmax  = 300;
t     = 3;
sigma = 0.01;

parfor jj=1:reps
    
    [FidCPMG(jj,:),FidXY2(jj,:)]=random_pi_pulse_errors(t,Nmax,sigma);
    
end

for N=1:Nmax
    
   meanF_CPMG(N)  = mean(FidCPMG(:,N));
   miner_CPMG(N)  = meanF_CPMG(N)-min(FidCPMG(:,N));
   maxer_CPMG(N)  = max(FidCPMG(:,N))-meanF_CPMG(N);
   
   meanF_XY2(N)  = mean(FidXY2(:,N));
   miner_XY2(N)  = meanF_XY2(N)-min(FidXY2(:,N));
   maxer_XY2(N)  = max(FidXY2(:,N))-meanF_XY2(N);
   
end

%% Plot as 2 subplots
close all

color=[24,116,205]/255;

color1 = [0,0,0];
color2 = [0,191,255]/255;%[152,245,255]/255;

subplot(2,1,1)

for jj=1:reps
   
    l=plot(1:Nmax,FidCPMG(jj,:),'color',[color1,0.1]);
    hold on
end

hold on

h=plot(1:Nmax,meanF_CPMG)
h.Marker='o';
h.MarkerFaceColor = color;
h.MarkerEdgeColor = color;
h.MarkerSize = 6;
h.LineWidth=2;
h.Color=color1;

xlabel('$N$','interpreter','latex')
ylabel('$\bar{F}$','interpreter','latex')
set(gca,'fontsize',24,'fontname','Microsoft Sans Serif')
title('CPMG')

set(gcf,'color','w')
xlim([1,Nmax])

set(gca, 'Box', 'on', 'TickDir', 'out', 'TickLength', [.003 .003], ...
    'XMinorTick', 'off', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
    'LineWidth', 0.5)


subplot(2,1,2)

for jj=1:reps
   
    l=plot(1:Nmax,FidXY2(jj,:),'color',[color1,0.1]);
    hold on
end
h=plot(1:Nmax,meanF_XY2)
h.Marker='o';
h.MarkerFaceColor = color;
h.MarkerEdgeColor = color;
h.MarkerSize = 6;
h.Color=color2;
h.LineWidth=2;

title('XY2')
xlabel('$N$','interpreter','latex')
ylabel('$\bar{F}$','interpreter','latex')
set(gca,'fontsize',24,'fontname','Microsoft Sans Serif')

set(gcf,'color','w')
xlim([1,Nmax])

set(gca, 'Box', 'on', 'TickDir', 'out', 'TickLength', [.003 .003], ...
    'XMinorTick', 'off', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
    'LineWidth', 0.5)





%% Plot with errorbars
close all



color1 = [0,0,0];
color2 = [0,191,255]/255;%[152,245,255]/255;

h=plot(1:Nmax,[meanF_CPMG;meanF_XY2])
h(1).Marker='s';
h(1).MarkerFaceColor = color1;
h(1).MarkerSize = 10;
h(2).Marker='o';
h(2).MarkerSize = 10;
h(2).MarkerFaceColor = color2;
h(1).LineWidth=2.5;
h(2).LineWidth=1.5;
h(1).Color=color1;
h(2).Color=color2;

hold on

for jj=1:reps
   
    l=plot(1:Nmax,FidCPMG(jj,:),'color',[color1,0.1]);
    
end

% hold on
% er1=errorbar(1:Nmax,h(1).YData,miner_CPMG,maxer_CPMG)
% er1.LineWidth=2.5;
% er1.Color=color1;
% er1.CapSize=8;
% 
% hold on
% er2=errorbar(1:Nmax,h(2).YData,miner_XY2,maxer_XY2)
% er2.LineWidth=1.5;
% er2.Color=color2;
% er2.CapSize=8;


xlabel('$N$','interpreter','latex')
ylabel('$F_{av}$','interpreter','latex')
legend('CPMG','XY2','location','best','color','none','edgecolor','none')
set(gca,'fontsize',24,'fontname','Microsoft Sans Serif')

set(gcf,'color','w')
xlim([1,Nmax])

set(gca, 'Box', 'on', 'TickDir', 'out', 'TickLength', [.003 .003], ...
    'XMinorTick', 'off', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
    'LineWidth', 0.5)


