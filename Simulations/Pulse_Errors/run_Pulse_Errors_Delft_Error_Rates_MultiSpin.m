function run_Pulse_Errors_Delft_Error_Rates_MultiSpin
%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
%
%Scirpt to add pulse errors for the Rx pulses based on the error rates
%reported in G. De Lange. “Quantum control and coherence of interacting spins in diamond”. PhD thesis. Delft University of Technology. (2012)
%
%This script calculates the effect of pulse errors on the M-tangling power
%of the target subspace for GHZ4 states, assuming the multi-spin scheme.
%Ref:arXiv:2302.05580
%--------------------------------------------------------------------------


close all

%Yellows={[255,185,15]/255,[255,215,0]/255,[238,238,0]/255};    
%Grays={[119,136,153]/255,[159,182,205]/255,[112,128,144]/255};
%Blues={[0,104,139]/255,[25,25,112]/255,[65,105,225]/255,[39,64,139]/255,[0,0,128]/255};
%Greens={[127,255,212]/255,[102,205,170]/255,[69,139,116]/255}  ;
%Reds={[205,55,0]/255,[205,38,38]/255,[139,0,0]/255};

LightBlue={[178,223,238]/255,[135,206,250]/255,[162,181,205]/255,[102,139,139]/255};

%============ Parameters ================================================

path        = '/Users/evatakou/Documents/MATLAB/Defect_Nuclear_GHZ_states/GHZ_states_Preparation/Simulations/GHZ_Data_of_MultiSpin_NEW/GHZ4_Multispin.mat';
load(path,"OUT");
[s0,s1,wL]  = load_fixed_params;
k           = 1;
cases       = length(OUT);
Nnuc        = length(OUT{1}.Target_Nuclei);
M           = Nnuc+1;
pref        = (2/3)^M;
pulse_error = -2/100;
%=========================================================================


epM_XY2_err    = zeros(1,length(cases));
epM_CPMG_ideal = zeros(1,length(cases));

for jj=1:cases
    
    At = OUT{jj}.A_Target;
    Bt = OUT{jj}.B_Target;
    t = OUT{jj}.Times;
    N = OUT{jj}.Iters;
    
    %=========== XY2 error ================================================
  
    temp = SuperClass_Sequences(wL,At,Bt,s0,s1,Nnuc,k,N);
    temp = temp.XY2_error_FULL(t,N,pulse_error);
    U    = temp.Uevol;
    
    temp            = SubClass_Ent_and_Fid;
    temp.Uval       = U;
    temp            = temp.MWayEp;
    epM_XY2_err(jj) = real(temp.epM_uni/pref);
    
    %=========== CPMG ideal ===============================================
  
    temp = SuperClass_Sequences(wL,At,Bt,s0,s1,Nnuc,k,N);
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

ylim([0.2,1])

legs={'CPMG ideal','XY2 with errors'};


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