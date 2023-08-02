function [p,minTau,tauPure,Tau_vec1,Tau_vec2,x_conv,y_conv]=Run_Minimal_3Tangle_Sequential

path = '/Users/evatakou/Documents/MATLAB/Defect_Nuclear_GHZ_states/GHZ_states_Preparation/Simulations/GHZ_Data_of_Sequential_NEW/GHZ3_Sequential.mat';
load(path,"OUT")

MaxCase =  length(OUT.Target_Nuclei);

parfor ii=1:MaxCase

    At    = OUT.A_Target{ii};
    Bt    = OUT.B_Target{ii};
    Aunw  = OUT.A_Unwanted{ii};
    Bunw  = OUT.B_Unwanted{ii};
    times = OUT.Opt_Unit_Times{ii};
    iters = OUT.Opt_Iters{ii};
    
    TEMP            = minimize_MTangle(At,Bt,Aunw,Bunw,times,iters);
    p(ii)           = TEMP.p;
    minTau(ii)      = TEMP.minTau;
    tauPure(ii)     = TEMP.PureTangle;
    Tau_vec1(ii)    = TEMP.tau_Vec1;
    Tau_vec2(ii)    = TEMP.tau_Vec2;
    
end

[p,indx]  = sort(p); 
Tau_vec1  = Tau_vec1(indx); %re-arrange based on sort order
Tau_vec2  = Tau_vec2(indx);
minTau    = minTau(indx);
tauPure   = tauPure(indx);
P         = zeros(length(p),2);

for ii=1:length(p)
    
   P(ii,1) = p(ii);
   P(ii,2) = minTau(ii);
   
end

[k,~]=convhull(P);

x_conv = P(k,1);
y_conv = P(k,2);




end



