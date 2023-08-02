function [Target_One_Tangle_Tol,Unwanted_One_Tangle_Tol,Infid_Tol,Time_Max,MaxRes,Time_Vals]=...
                     Tolerances_Multispin(GHZsize)

MaxRes=10; 

if GHZsize==3
    
    Time_Max                = 2000;
    Target_One_Tangle_Tol   = 0.9;
    Unwanted_One_Tangle_Tol = 0.1;
    Infid_Tol               = 0.11;
    Time_Vals               = 500;  
    
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