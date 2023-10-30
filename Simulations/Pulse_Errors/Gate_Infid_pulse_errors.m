function Gate_Infid_pulse_errors(Option)

%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 29, 2023
%--------------------------------------------------------------------------


close all;


%------------- Definitions ------------------------------------------------
s0=0;
s1=-1;

X=[0 1; 1 0];
Y=[0 -1i ; 1i 0];
Z=[1 0 ; 0 -1];
Ze=[s0 0 ; 0 s1];
Rx=@ (err) expm(-1i*(pi+err)/2*X);
Ry=@ (err) expm(-1i*(pi+err)/2*Y);


ket0 = [1;0];
ket1 = [0;1];
ket2 = [1;1]/sqrt(2);
ket3 = [1;-1]/sqrt(2);
ket4 = [1;1i]/sqrt(2);
ket5 = [1;-1i]/sqrt(2);

kets = {ket0,ket1,ket2,ket3,ket4,ket5};

for jj=1:length(kets)

    rho0s{jj}=kets{jj}*kets{jj}';

end

%--------------------------------------------------------------------------
c         = 2*pi*1e-3;       %Conversion factor to angular frequency (in kHz)  
B0        = 403;             %External B-field in Gauss
B0        = B0*1e-4;         %External B-field in Tesla 
t         = 3;               %\mu s
gamma_e   = 1.760859*10^11 ; %rad/T*s

we    = gamma_e * B0 ;          %Larmor freq of electron rad/s
we    = we*1e-6;                %rad*MHz                     
we    = we*1e3;                 %rad*KHz
we    = we/(2*pi);              %kHz       
we    = c*we;                   %kHz*rad

Hfree = we*Ze;     
Uf    = @(t) expm(-1i*t*Hfree); 
%--------------------------------------------------------------------------

UCPMG=@(t,err) Uf(t/4)*Rx(err(1))*Uf(t/2)*Rx(err(2))*Uf(t/4); %CPMG unit
UXY2=@(t,err)  Uf(t/4)*Ry(err(1))*Uf(t/2)*Rx(err(2))*Uf(t/4); %XY2 unit

%--------------------------------------------------------------------------

Nmax  = 300;
reps  = 1000;
sigma = 0.01;

switch Option
    
    case 'Systematic'
        
        for N=1:Nmax
            
            %U0_CPMG = UCPMG(t,[0,0])^N;
            U0_XY2  = UXY2(t,[0,0])^N;
            
            %Uerr_CPMG_0_02 = UCPMG(t,[0.02,0.02])^N;
            Uerr_XY2_0_02  = UXY2(t,[2/100,2/100])^N;
            Uerr_XY2_0_08  = UXY2(t,[8/100,8/100])^N;
            
            
            %F_CPMG(N)     = calculate_fid(rho0s,U0_CPMG,Uerr_CPMG_0_02);
            F_XY2_0_02(N) = calculate_fid(rho0s,U0_XY2,Uerr_XY2_0_02);
            F_XY2_0_08(N) = calculate_fid(rho0s,U0_XY2,Uerr_XY2_0_08);
            
            
        end
        
        h=semilogy(1:Nmax,[1-F_XY2_0_02;1-F_XY2_0_08],'linewidth',2); 
        
        
        h(1).Marker='o';

        legend('XY2 $\epsilon=0.02$','XY2 $\epsilon=0.08$',...
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
        
    case 'Random'
        
        parfor ii=1:reps
        
            [FidCPMG(ii,:),FidXY2(ii,:)]=one_realization_random_errors(UCPMG,UXY2,Nmax,t,sigma,rho0s);
        
        end
        
        %Get the average, and deviation from mean
        
        for N=1:Nmax
            
           meanF_CPMG(N)  = mean(FidCPMG(:,N));
           miner_CPMG(N)  = meanF_CPMG(N)-min(FidCPMG(:,N));
           maxer_CPMG(N)  = max(FidCPMG(:,N))-meanF_CPMG(N);

           meanF_XY2(N)  = mean(FidXY2(:,N));
           miner_XY2(N)  = meanF_XY2(N)-min(FidXY2(:,N));
           maxer_XY2(N)  = max(FidXY2(:,N))-meanF_XY2(N);
            
        end
        
        %-------- Plot the results ---------------------------------------
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
        
        
        
end


end


function F=calculate_fid(rho0s,U0,Uerr)


for jj=1:length(rho0s)
    
    rho0 = rho0s{jj};
    
    rho_err = Uerr*rho0*Uerr';
    
    f(jj)=real( trace( U0*rho0*U0'*rho_err ) );
    
    
end

F=sum(f)/6;



end

function [FidCPMG,FidXY2]=one_realization_random_errors(UCPMG,UXY2,Nmax,t,sigma,rho0s)

Uer_CPMG = eye(2);
Uer_XY2  = eye(2);


for N=1:Nmax

    U0_CPMG  = UCPMG(t,[0,0])^N;
    U0_XY2   = UXY2(t,[0,0])^N;

    r        = normrnd(0,sigma,[1,2]);
    
    Uer_CPMG = Uer_CPMG*UCPMG(t,r);
    Uer_XY2  = Uer_XY2*UXY2(t,r);
    

    FidCPMG(N) = calculate_fid(rho0s,U0_CPMG,Uer_CPMG);
    FidXY2(N)  = calculate_fid(rho0s,U0_XY2,Uer_XY2);

end


end


