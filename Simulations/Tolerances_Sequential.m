function [Target_One_Tangle_Tol,Unwanted_One_Tangle_Tol,Infid_Tol,Time_Max,Individual_Time_Max,MaxRes]=...
                     Tolerances_Sequential(GHZsize)
%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
%
%Script to set the tolerances for the optimization of the sequential 
%entangling scheme optimization to reproduce the results
%of Ref. https://arxiv.org/abs/2302.05580.
%--------------------------------------------------------------------------
%
%Input: GHZsize>=3, which is the size of the GHZ state.
%--------------------------------------------------------------------------

MaxRes=18;  %Max resonance index

if GHZsize==3
    
    Time_Max                = 2000; %Total time of composite sequence (\mu s)
    Individual_Time_Max     = 2000; %Time of individual sequence (\mu s)
    Unwanted_One_Tangle_Tol = 0.1;  %Tolerance of unwanted one-tangles 
    Target_One_Tangle_Tol   = 0.99; %Tolerance of target one-tangles
    Infid_Tol               = 0.1;  %Tolerance for infidelity of composite gate
    
elseif GHZsize==4
    
    Time_Max                = 2000; 
    Individual_Time_Max     = 2000;
    Unwanted_One_Tangle_Tol = 0.1;
    Target_One_Tangle_Tol   = 0.99;
    Infid_Tol               = 0.1;

elseif GHZsize==5
    
    Time_Max                = 2300; 
    Individual_Time_Max     = 2300;
    Unwanted_One_Tangle_Tol = 0.1;
    Target_One_Tangle_Tol   = 0.9;
    Infid_Tol               = 0.1;    
    
elseif GHZsize==6

    Time_Max                = 2500; 
    Individual_Time_Max     = 2500;
    Unwanted_One_Tangle_Tol = 0.12;
    Target_One_Tangle_Tol   = 0.9;
    Infid_Tol               = 0.11;    
    
elseif GHZsize==7

    Time_Max                = 3300; 
    Individual_Time_Max     = 3300;
    Unwanted_One_Tangle_Tol = 0.12;
    Target_One_Tangle_Tol   = 0.9;
    Infid_Tol               = 0.12;    
    
elseif GHZsize==8
    
    Time_Max                = 3700; 
    Individual_Time_Max     = 2000;
    Unwanted_One_Tangle_Tol = 0.12;
    Target_One_Tangle_Tol   = 0.9;
    Infid_Tol               = 0.13;    

    
elseif GHZsize==9
    
    Time_Max                = 4000; 
    Individual_Time_Max     = 1400;
    Unwanted_One_Tangle_Tol = 0.15;
    Target_One_Tangle_Tol   = 0.85;
    Infid_Tol               = 0.13;     
    
elseif GHZsize==10
    
    Time_Max                = 4000; 
    Individual_Time_Max     = 1400;
    Unwanted_One_Tangle_Tol = 0.22;
    Target_One_Tangle_Tol   = 0.87;
    Infid_Tol               = 0.19;
    
end


end