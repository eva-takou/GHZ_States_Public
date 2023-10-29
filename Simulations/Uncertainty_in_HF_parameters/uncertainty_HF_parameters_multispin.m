function uncertainty_HF_parameters_multispin(p)
%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 29, 2023
%--------------------------------------------------------------------------
%
%Script to add uncertainty in HF parameters by shifting the HF components
%of target nuclei by constant p value. This script compares the target subspace M-tangling
%power obtained by no uncertainty in HF parameters, with the one with
%certain HF parameters. It tests the multipsin scheme, see
%arXiv:2302.05580.
%
%Input: p amount to shift the HF parameters (units assumed to be kHz)
%--------------------------------------------------------------------------



close all;

GHZsizeMax=6;

%============ Parameters ================================================

[s0,s1,wL]=load_fixed_params;

%=========================================================================
tiledlayout('flow')

for ll=3:GHZsizeMax
    
    GHZsize=ll;

    if GHZsize==3
        load('/Users/evatakou/Documents/MATLAB/Defect_Nuclear_GHZ_states/GHZ_states_Preparation/Simulations/GHZ_Data_of_MultiSpin_NEW/GHZ3_Multispin.mat',"OUT")
    elseif GHZsize==4
        load('/Users/evatakou/Documents/MATLAB/Defect_Nuclear_GHZ_states/GHZ_states_Preparation/Simulations/GHZ_Data_of_MultiSpin_NEW/GHZ4_Multispin.mat',"OUT")
    elseif GHZsize==5
        load('/Users/evatakou/Documents/MATLAB/Defect_Nuclear_GHZ_states/GHZ_states_Preparation/Simulations/GHZ_Data_of_MultiSpin_NEW/GHZ5_Multispin.mat',"OUT")
    elseif GHZsize==6
        load('/Users/evatakou/Documents/MATLAB/Defect_Nuclear_GHZ_states/GHZ_states_Preparation/Simulations/GHZ_Data_of_MultiSpin_NEW/GHZ6_Multispin.mat',"OUT")
    end
    
    
    cases = length(OUT);

    Nnuc = 1;
    k    = 1;
    epM  = zeros(1,cases);
    ep_nuc  = zeros(cases,length(OUT{1}.A_Target));
    epM_uncertain = zeros(1,cases);
    for jj=1:cases

        At = OUT{jj}.A_Target;
        Bt = OUT{jj}.B_Target;
        t  = OUT{jj}.Times;
        N  = OUT{jj}.Iters;

       %For the sequential we need to compose the CPMG sequences.

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
    
      clear phi0
      clear phi1
      clear n0
      clear n1
      clear ep_nuc      

      At = OUT{jj}.A_Target;
      Bt = OUT{jj}.B_Target;

      %At = At*(1+p);
      %Bt = Bt*(1+p);
      At = At+p;
      Bt = Bt+p;

      for indx=1:length(At)

          temp = SubClass_U4Operations(wL,At(indx),Bt(indx),s0,s1,Nnuc,k,1);
          temp = temp.CPMG(t,N);
          temp = temp.Rot_Angles_And_Axes;

          n0{indx}   = temp.axes{1};
          n1{indx}   = temp.axes{2};
          phi0(indx) = temp.angles{1};
          phi1(indx) = temp.angles{2};

          n0n1 = dot(n0{indx},n1{indx});
          ep_nuc(jj,indx) = 1-(cos(phi0(indx)/2)*cos(phi1(indx)/2)+n0n1*sin(phi0(indx)/2)*sin(phi1(indx)/2))^2;

      end   

      epM_uncertain(jj) = prod(ep_nuc(jj,:));

    end

    
    nexttile
    h=bar((epM-epM_uncertain));
    h(1).FaceColor=[0,104,139]/255;
    
    if length(epM)==15
       set(gca,'xtick',1:3:15) 
    end
    
    clear phi0
    clear phi1
    clear n0
    clear n1
    clear onetangle
    
    clear epM
    clear epM_uncertain


    xlabel('Case #')
    ylabel('$\Delta \epsilon_{p,M}(U)$','interpreter','latex')
    set(gcf,'color','w')
    set(gca,'fontsize',18,'fontname','Microsoft Sans Serif')
    title(strcat('GHZ_',num2str(ll)))
    set(gca, 'Box', 'on', 'TickDir', 'out', 'TickLength', [.003 .003], ...
    'XMinorTick', 'off', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
    'LineWidth', 0.5)   
    
end





end