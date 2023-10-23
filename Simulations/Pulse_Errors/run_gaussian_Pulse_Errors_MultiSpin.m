function run_gaussian_Pulse_Errors_MultiSpin(variance)
%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
%
%Script to add random over-/under-rotation errors for the Rx(pi) pulses,
%sampled from a normal distribution and calculate the M-tangling power of
%the target subspace. Considers the multispin protocol of arXiv:2302.05580,
%and GHZ4 states.
%
%Input: variance of normal distribution
%
%--------------------------------------------------------------------------
close all

%Yellows={[255,185,15]/255,[255,215,0]/255,[238,238,0]/255};    
%Grays={[119,136,153]/255,[159,182,205]/255,[112,128,144]/255};
%Blues={[0,104,139]/255,[25,25,112]/255,[65,105,225]/255,[39,64,139]/255,[0,0,128]/255};
%Greens={[127,255,212]/255,[102,205,170]/255,[69,139,116]/255}  ;
%Reds={[205,55,0]/255,[205,38,38]/255,[139,0,0]/255};

LightBlue={[178,223,238]/255,[135,206,250]/255,[162,181,205]/255,[102,139,139]/255};

%============ Parameters ================================================

[s0,s1,wL]  = load_fixed_params;
path        = '/Users/evatakou/Documents/MATLAB/Defect_Nuclear_GHZ_states/GHZ_states_Preparation/Simulations/GHZ_Data_of_MultiSpin_NEW/GHZ4_Multispin.mat';
load(path,"OUT");
k           = 1;
cases       = length(OUT);
Nnuc        = length(OUT{1}.Target_Nuclei);
M           = Nnuc+1;
pref        = (2/3)^M;

%=========================================================================

MaxRep = 500; %~20mins

D=parallel.pool.Constant(OUT);

tic
parfor ll=1:MaxRep
    
   for jj=1:cases
       
       At = D.Value{jj}.A_Target;
       Bt = D.Value{jj}.B_Target;
   
       t = D.Value{jj}.Times;
       N = D.Value{jj}.Iters;

       %========= CPMG error ==============================================
       
       temp = SuperClass_Sequences(wL,At,Bt,s0,s1,Nnuc,k,N);
       temp = temp.CPMG_error_variance(t,N,variance);
       U    = temp.Uevol;
       
       temp                     = SubClass_Ent_and_Fid;
       temp.Uval                = U;
       temp                     = temp.MWayEp;
       epM_temp_CPMG_err(jj,ll) = real(temp.epM_uni/pref)
       
       %========= XY2 error ===============================================
       
       temp = SuperClass_Sequences(wL,At,Bt,s0,s1,Nnuc,k,N);
       temp = temp.XY2_error_variance(t,N,variance);
       U    = temp.Uevol;
       
       temp                    = SubClass_Ent_and_Fid;
       temp.Uval               = U;
       temp                    = temp.MWayEp;
       epM_temp_XY2_err(jj,ll) = real(temp.epM_uni/pref)
        
   end
    
end
toc

epM_CPMG_ideal   = zeros(1,cases);
epM_CPMG_err     = zeros(1,cases);
err_bar_min_CPMG = zeros(1,cases);
err_bar_max_CPMG = zeros(1,cases);
epM_XY2_err      = zeros(1,cases);
err_bar_min_XY2  = zeros(1,cases);
err_bar_max_XY2  = zeros(1,cases);

for jj=1:cases

       At = D.Value{jj}.A_Target;
       Bt = D.Value{jj}.B_Target;
   
       t = D.Value{jj}.Times;
       N = D.Value{jj}.Iters;
    
       temp = SuperClass_Sequences(wL,At,Bt,s0,s1,Nnuc,k,N);
       temp = temp.CPMG_error_variance(t,N,0*variance);
       U    = temp.Uevol;
       
       temp                    = SubClass_Ent_and_Fid;
       temp.Uval               = U;
       temp                    = temp.MWayEp;
       epM_CPMG_ideal(jj)   = real(temp.epM_uni/pref);
    
       epM_CPMG_err(jj)     = mean(epM_temp_CPMG_err(jj,:));
       err_bar_min_CPMG(jj) = abs(min(epM_temp_CPMG_err(jj,:))-epM_CPMG_err(jj));
       err_bar_max_CPMG(jj) = abs(max(epM_temp_CPMG_err(jj,:))-epM_CPMG_err(jj));
    
       epM_XY2_err(jj)     = mean(epM_temp_XY2_err(jj,:));
       err_bar_min_XY2(jj) = abs(epM_XY2_err(jj)-min(epM_temp_XY2_err(jj,:)));
       err_bar_max_XY2(jj) = abs(max(epM_temp_XY2_err(jj,:))-epM_XY2_err(jj));

end


model_series = 1:cases;
err_min      = [err_bar_min_XY2;err_bar_min_CPMG];
err_max      = [err_bar_max_XY2;err_bar_max_CPMG];
h=bar(model_series,[epM_CPMG_ideal;epM_XY2_err;epM_CPMG_err]);

h(1).FaceColor=[0,104,139]/255;
h(2).FaceColor=LightBlue{1};
h(3).FaceColor=[0.89,0.63,0.31];

hold on
% Find the number of groups and the number of bars in each group
ngroups = length(1:cases);
nbars   = 3;

% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 2:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*(i)-1) * groupwidth / (2*nbars);
    er=errorbar(x,h(i).YData, err_min(i-1,:),err_max(i-1,:), 'k', 'linestyle', 'none');
    er.LineStyle='none';
    er.Color=[0,0,0];
    er.LineWidth=1.5;
    er.CapSize=10;
    
end
hold off





ylim([0.5,1])

legs={'CPMG ideal',strcat('XY2,~$\sigma=$',num2str(variance)),...
        strcat('CPMG,~$\sigma=$',num2str(variance))};


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