function [p,minTau,tauPure,Tau_vec1,Tau_vec2,x_conv,y_conv]=Run_Minimal_3Tangle_Sequential
%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
%
%Function for convex hull construction of the 3-tangle of a mixed state (it has rank 2).
%Uses methods initially introduced in
%https://dx.doi.org/10.1103/PhysRevA.77.032310, and is used to produce the
%results in arXiv:2302.05580.
%
%
%Output: p: largest e-val of density matrix
%        minTau: the minimal value of M-tangle of the mixed state
%        tauPure: the M-tangle assuming a pure state (w/o tracing out the unwanted subsystems)    
%        Tau_vec1: the M-tangle of 1st component of mixed state 
%        Tau_vec2: the M-tangle of 2nd component of mixed state
%        x_conv/y_conv: the convex hull of the mixed state tangle as a
%        function of largest eigenvalue.


path = '/Users/evatakou/Documents/MATLAB/GHZ_States_Public/Simulations/GHZ_Data_of_Sequential_FINAL/GHZ3_Sequential.mat';
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



