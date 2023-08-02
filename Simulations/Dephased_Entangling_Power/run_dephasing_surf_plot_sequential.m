function run_dephasing_surf_plot_sequential

%============ Parameters ================================================
[s0,s1,wL]   = load_fixed_params;
path         = '/Users/evatakou/Documents/MATLAB/Defect_Nuclear_GHZ_states/GHZ_states_Preparation/Simulations/GHZ_Data_of_Sequential_NEW/GHZ3_Sequential.mat';
load(path,"OUT")
Nnuc         = 1;
k            = 1;
lambda0      = 0.8:0.001:1;
lambda1      = 1-lambda0;
theta_coeffs = 1/500:0.0001:1/50; %1/mus
caseNum      = 1;
d            = 2;
M            = length(OUT.A_Target{caseNum})+1;
pref         = (d/(d+1))^M;
%=========================================================================

At       = OUT.A_Target{caseNum};
Bt       = OUT.B_Target{caseNum};
t        = OUT.Opt_Unit_Times{caseNum};
N        = OUT.Opt_Iters{caseNum};
n0       = cell(1,length(At));
n1       = n0;
phi0     = zeros(1,length(At));
phi1     = phi0;
epM_Deph = zeros(length(lambda0),length(theta_coeffs));
   
for indx=1:length(At)

  test = SubClass_U4Operations(wL,At(indx),Bt(indx),s0,s1,Nnuc,k,1);
  test = test.CPMG(t,N);
  test = test.Rot_Angles_And_Axes;

  n0{indx}   = test.axes{1};
  n1{indx}   = test.axes{2};
  phi0(indx) = test.angles{1};
  phi1(indx) = test.angles{2};

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