function run_systematic_Pulse_Errors_Sequential(pulse_error)
%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 28, 2023
%--------------------------------------------------------------------------
%
%Scirpt to add pulse errors for the Rx pulses and find the effect on the M-tangling power
%of GHZ4 states, using the sequential scheme. Ref:arXiv:2302.05580
%--------------------------------------------------------------------------
%
%Input: Pulse error is systematic pulse error of Rx(pi) (% of pi)
%--------------------------------------------------------------------------

close all
%Yellows={[255,185,15]/255,[255,215,0]/255,[238,238,0]/255};    
%Grays={[119,136,153]/255,[159,182,205]/255,[112,128,144]/255};
%Blues={[0,104,139]/255,[25,25,112]/255,[65,105,225]/255,[39,64,139]/255,[0,0,128]/255};
%Greens={[127,255,212]/255,[102,205,170]/255,[69,139,116]/255}  ;
%Reds={[205,55,0]/255,[205,38,38]/255,[139,0,0]/255};

LightBlue={[178,223,238]/255,[135,206,250]/255,[162,181,205]/255,[102,139,139]/255};

%============ Parameters ================================================
path        = '/Users/evatakou/Documents/MATLAB/GHZ_States_Public/Simulations/GHZ_Data_of_Sequential_FINAL/GHZ4_Sequential.mat';
load(path,"OUT");
[s0,s1,wL]  = load_fixed_params;
k           = 1;
cases       = length(OUT.Target_Nuclei);
Nnuc        = length(OUT.Target_Nuclei{1});
M           = Nnuc+1;
pref        = (2/3)^M;
%=========================================================================

epM_CPMG_err   = zeros(1,length(cases));
epM_XY2_err    = epM_CPMG_err;
epM_CPMG_ideal = epM_CPMG_err;

for jj=1:cases
    
    At = OUT.A_Target{jj};
    Bt = OUT.B_Target{jj};
    t  = OUT.Opt_Unit_Times{jj};
    N  = OUT.Opt_Iters{jj};
   
    %=================== CPMG with pulse errors ===========================
    temp = SuperClass_Sequences(wL,At,Bt,s0,s1,Nnuc,k,1);  
    temp = temp.CPMG_error(t,N,pulse_error);  
    U    = temp.Uevol;
  
    temp             = SubClass_Ent_and_Fid;
    temp.Uval        = U;
    temp             = temp.MWayEp;
    epM_CPMG_err(jj) = real(temp.epM_uni/pref);

    %=================== XY2 with pulse errors ============================
    temp = SuperClass_Sequences(wL,At,Bt,s0,s1,Nnuc,k,1);
    temp = temp.XY2_error(t,N,pulse_error);
    U    = temp.Uevol;  
    
    temp             = SubClass_Ent_and_Fid;
    temp.Uval        = U;
    temp             = temp.MWayEp;
    epM_XY2_err(jj)  = real(temp.epM_uni/pref);
  
    %============== CPMG ideal ===========================================
    temp = SuperClass_Sequences(wL,At,Bt,s0,s1,Nnuc,k,1);
    temp = temp.CPMG_error(t,N,0*pulse_error);  
    U    = temp.Uevol;
    
    temp               = SubClass_Ent_and_Fid;
    temp.Uval          = U;
    temp               = temp.MWayEp;
    epM_CPMG_ideal(jj) = real(temp.epM_uni/pref);
    
   
end

bar(epM_CPMG_ideal,'Facecolor',[0,104,139]/255,'FaceAlpha',1)
hold on
bar(epM_XY2_err,'Facecolor',LightBlue{1},'FaceAlpha',1,'BarWidth',0.4)
hold on
bar(epM_CPMG_err,'Facecolor',[0.89,0.63,0.31],'FaceAlpha',1,'BarWidth',0.3) %Greens 2

ylim([0.5,1])

legs={'CPMG ideal',strcat('XY2,~',num2str(pulse_error*100),'\%',' error'),...
        strcat('CPMG,~',num2str(pulse_error*100),'\%',' error')};


legend(legs,'color','none',...
    'edgecolor','none','location','best','interpreter','latex','NumColumns',5)

set(gca, 'Box', 'on', 'TickDir', 'out', 'TickLength', [.003 .003], ...
    'XMinorTick', 'off', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
    'LineWidth', 0.5)

xlabel('Case #')
ylabel('$\epsilon_{p,M}(U)$','interpreter','latex')
set(gcf,'color','w')
set(gca,'fontsize',24,'fontname','Microsoft Sans Serif')


end


