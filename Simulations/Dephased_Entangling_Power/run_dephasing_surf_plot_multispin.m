function run_dephasing_surf_plot_multispin
%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 28, 2023
%--------------------------------------------------------------------------
%Function to plot the M-tangling power when the electron undergoes
%dephasing [see Ref. arXiv:2302.05580], and for the multispin protocol.
%
%The dephasing channel assumes the Kraus operators: 
%
%K0 = sqrt(lambda0) [exp(i*theta*t) 0]
%                   [0   exp(-i*theta*t)]
%
%K1 = sqrt(lambda1) [exp(-i*theta*t) 0]
%                   [0   exp(i*theta*t)]
%
%1/theta [in mus] describes time-scale of dephasing. It holds that
%lambda0+lambda1=1. [Kraus operators from Ref. https://dx.doi.org/10.1103/PhysRevA.106.022414]

%============ Parameters ================================================

[s0,s1,wL]  = load_fixed_params;
path        = '/Users/evatakou/Documents/MATLAB/GHZ_States_Public/Simulations/GHZ_Data_of_MultiSpin_NEW/GHZ3_Multispin.mat';
load(path,"OUT")
Nnuc         = 1;
k            = 1;
lambda0      = 0.8:0.001:1;
lambda1      = 1-lambda0;
theta_coeffs = 1/500:0.0001:1/50; %1/mus
caseNum      = 1;
d            = 2;
M            = length(OUT{caseNum}.A_Target)+1;
pref         = (d/(d+1))^M;
%=========================================================================
 
At       = OUT{caseNum}.A_Target;
Bt       = OUT{caseNum}.B_Target;
t        = OUT{caseNum}.Times;
N        = OUT{caseNum}.Iters;
n0       = cell(1,length(At));
n1       = n0;
phi0     = zeros(1,length(At));
phi1     = zeros(1,length(At));
epM_Deph = zeros(length(lambda0),length(theta_coeffs));

for indx=1:length(At)

  temp = SubClass_U4Operations(wL,At(indx),Bt(indx),s0,s1,Nnuc,k,1);
  temp = temp.CPMG(t,N);
  temp = temp.Rot_Angles_And_Axes;

  n0{indx}   = temp.axes{1};
  n1{indx}   = temp.axes{2};
  phi0(indx) = temp.angles{1};
  phi1(indx) = temp.angles{2};

end


for kk=1:length(lambda0)

  for ll=1:length(theta_coeffs)

      temp = SubClass_Ent_and_Fid;
      temp = temp.dephased_epM_CR(n0,n1,phi0,phi1,lambda0(kk),lambda1(kk),theta_coeffs(ll),t,N);
      epM_Deph(kk,ll) = temp.epM_dephased;
      
  end

end


h=surf(lambda0,theta_coeffs.^(-1),epM_Deph'/pref);

h.EdgeColor='none';

xlabel('$\lambda_0$','interpreter','latex')
ylabel('$\theta^{-1}$ ($\mu$s)','interpreter','latex')
set(gcf,'color','w')
set(gca,'fontsize',24,'fontname','Microsoft Sans Serif')
view(0,90)

shading interp
colorbar
colormap(brewermap([],'PuBu'))
set(gca, 'Box', 'on', 'TickDir', 'out', 'TickLength', [.003 .003], ...
    'XMinorTick', 'off', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
    'LineWidth', 0.5)

ylim([min(1./theta_coeffs),max(1./theta_coeffs)])


end