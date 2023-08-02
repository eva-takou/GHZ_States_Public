function run_dephasing_bar_plot_sequential(lambda0)

%========== Colors ========================================================

%Blues={[0,104,139]/255,[25,25,112]/255,[65,105,225]/255,[39,64,139]/255,[0,0,128]/255};

%Magentas={[0.4940 0.1840 0.5560],[147,112,219]/255};

%Purples={[204,153,255]/255,[153,0,153]/255,[102,0,51]/255,...
%        [51,0,102]/255};
    
Yellows={[255,185,15]/255,[255,215,0]/255,[238,238,0]/255};    

Greens={[127,255,212]/255,[102,205,170]/255,[69,139,116]/255}  ;

%Grays={[119,136,153]/255,[159,182,205]/255,[112,128,144]/255};

Reds={[205,55,0]/255,[205,38,38]/255,[139,0,0]/255};

%OtherColors={[205,179,139]/255};

colorComb={Yellows{2},Greens{2},Reds{1}};

%============ Parameters ================================================
[s0,s1,wL]  = load_fixed_params;
theta_coeff = [1/400,1/100,1/50]; %1/mus
lambda1     = 1-lambda0;
path        = '/Users/evatakou/Documents/MATLAB/Defect_Nuclear_GHZ_states/GHZ_states_Preparation/Simulations/GHZ_Data_of_Sequential_NEW/GHZ3_Sequential.mat';
load(path,"OUT")
cases       = length(OUT.Target_Nuclei);
Nnuc        = 1;
k           = 1;
d           = 2;
M           = length(OUT.A_Target{1})+1;
pref        = (d/(d+1))^M;
%=========================================================================

epM      = zeros(1,cases);
epM_Deph = zeros(cases,length(theta_coeff));
ep_nuc   = zeros(cases,length(OUT.A_Target{1}));
n0       = cell(length(OUT.A_Target{1}));
n1       = cell(length(OUT.A_Target{1}));
phi0     = zeros(1,length(OUT.A_Target{1}));
phi1     = zeros(1,length(OUT.A_Target{1}));

for jj=1:cases
    
    At = OUT.A_Target{jj};
    Bt = OUT.B_Target{jj};
    t  = OUT.Opt_Unit_Times{jj};
    N  = OUT.Opt_Iters{jj};
   
    for indx=1:length(At)
      
        temp = SubClass_U4Operations(wL,At(indx),Bt(indx),s0,s1,Nnuc,k,1);
        temp = temp.CPMG(t,N);
        temp = temp.Rot_Angles_And_Axes;
      
        n0{indx}   = temp.axes{1};
        n1{indx}   = temp.axes{2};
        phi0(indx) = temp.angles{1};
        phi1(indx) = temp.angles{2};
        n0n1       = dot(n0{indx},n1{indx});
      
        ep_nuc(jj,indx) = 1-(cos(phi0(indx)/2)*cos(phi1(indx)/2)+n0n1*sin(phi0(indx)/2)*sin(phi1(indx)/2))^2;
  
    end
    
    epM(jj)      = prod(ep_nuc(jj,:));
  
    for ll=1:length(theta_coeff)
  
        temp  = SubClass_Ent_and_Fid;
        temp  = temp.dephased_epM_CR(n0,n1,phi0,phi1,lambda0,lambda1,theta_coeff(ll),t,N);
        epM_Deph(jj,ll) = temp.epM_dephased;
        
    end
  
end

subplot(2,1,1)

h=bar([epM;epM_Deph'./pref]');
h(1).FaceColor=[0,104,139]/255;
h(1).EdgeColor='k';
legs = {'$\epsilon_{p,M}(U)$'};

cnt=0;
for ll=1:length(theta_coeff)

    if ll<4
        cnt=cnt+1;
        h(ll+1).FaceColor=colorComb{cnt};
        h(ll+1).EdgeColor='k';
    end
    legs=[legs,strcat('$\epsilon_{p,M}(\mathcal{E}_{deph}), \theta^{-1}=$',num2str(1/theta_coeff(ll)),'$\mu$s')];
    
end

xlabel('Case #')
ylabel('$\epsilon_{p,M}$','interpreter','latex')
set(gcf,'color','w')

legend(legs,'color','none','edgecolor','none',...
            'location','best','interpreter','latex','NumColumns',4)

set(gca, 'Box', 'on', 'TickDir', 'out', 'TickLength', [.003 .003], ...
         'XMinorTick', 'off', 'YMinorTick', 'on', 'YGrid', 'on', ...
         'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
         'LineWidth', 0.5)

set(gca,'fontsize',24,'fontname','Microsoft Sans Serif')
ylim([0.5,1])
set(gca,'ytick',[0.6,0.8,1])
end