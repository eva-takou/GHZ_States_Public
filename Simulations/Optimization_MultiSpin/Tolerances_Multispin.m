function [Target_One_Tangle_Tol,Unwanted_One_Tangle_Tol,Infid_Tol,Time_Max,MaxRes,Time_Vals]=...
                     Tolerances_Multispin(GHZsize)
%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
%
%Script to set the tolerances for the optimization of the multispin 
%entangling scheme optimization to reproduce the results
%of Ref. https://arxiv.org/abs/2302.05580.
%--------------------------------------------------------------------------
%
%Input: GHZsize>=3, which is the size of the GHZ state.
%--------------------------------------------------------------------------

MaxRes=10; %Max resonance index

if GHZsize==3
    
    Time_Max                = 2000; %Total time of sequence (\mu s)
    Target_One_Tangle_Tol   = 0.9;  %Tolerance of target one-tangles
    Unwanted_One_Tangle_Tol = 0.1;  %Tolerance of unwanted one-tangles 
    Infid_Tol               = 0.11; %Tolerance for infidelity of multispin gate
    Time_Vals               = 500;  %# of time values to sample (determines time-step)
    
elseif GHZsize==4
    
    Time_Max                = 1800;
    Target_One_Tangle_Tol   = 0.9;
    Unwanted_One_Tangle_Tol = 0.12;
    Infid_Tol               = 0.101;
    Time_Vals               = 500;
    
elseif GHZsize==5
    
    Time_Max                = 2300;
    Target_One_Tangle_Tol   = 0.84;
    Unwanted_One_Tangle_Tol = 0.13;
    Infid_Tol               = 0.11;
    Time_Vals               = 400;  

elseif GHZsize==6

    Time_Max                = 2500;
    Target_One_Tangle_Tol   = 0.88;
    Unwanted_One_Tangle_Tol = 0.12;
    Infid_Tol               = 0.13;
    Time_Vals               = 200;
    
elseif GHZsize==7
    
    Time_Max                = 2800;
    Target_One_Tangle_Tol   = 0.85;
    Unwanted_One_Tangle_Tol = 0.15;
    Infid_Tol               = 0.13;
    Time_Vals               = 200;
    
elseif GHZsize==8
    
    Time_Max                = 3000;
    Target_One_Tangle_Tol   = 0.85;
    Unwanted_One_Tangle_Tol = 0.15;
    Infid_Tol               = 0.15;
    Time_Vals               = 200;
    
elseif GHZsize==9
    
    Time_Max                = 3000;
    Target_One_Tangle_Tol   = 0.82;
    Unwanted_One_Tangle_Tol = 0.15;
    Infid_Tol               = 0.15;
    Time_Vals               = 200;
    
end





end