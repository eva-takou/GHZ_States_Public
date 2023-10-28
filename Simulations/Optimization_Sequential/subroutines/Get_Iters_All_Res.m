function OUT=Get_Iters_All_Res(Spin_Indx,Res_Max,wL,A,B,s0,s1,Number_Of_Maxima,Time_Max,Unwanted_One_Tangle_Tol)
%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 27, 2023
%--------------------------------------------------------------------------
%
%Output all possible Iters and optimal times for up to a max resonance and
%for a specific nuclear spin labeled by Spin_Indx.
%Input: Spin_Indx \in [1,length(A)]
%       Res_Max: Maximum resonance to search for
%       wL,A,B: Larmor frequency and HF parameters of 27 nuclei
%       s0,s1: Electron's projections
%       Number_Of_Maxima: how many maxima of one-tangle to keep 
%       Time_Max: Time constraint
%       Unwanted_One_Tangle_Tol: Tolerance of unwanted one-tangles.
%Output: Optimal unit times and iterations.


cnt=0;
CNT=0;

for Res=1:Res_Max
    
   [t,N]=Get_Iters(Spin_Indx,Res,wL,A,B,s0,s1,Number_Of_Maxima,Time_Max,Unwanted_One_Tangle_Tol) ;
    
    if ~isnan(t) %This resonance satisfied one-tangle constraints we imposed
        
       for ll=1:length(N)  %For each resonance there are multiple iterations we keep
       
           cnt            = cnt+1; 
           OUT.Times(cnt) = t;
           OUT.Iters(cnt) = N(ll);
       
       end
       
    else
        
        CNT=CNT+1;
       
    end
                   
    
end

if CNT==Res_Max %We didnt satisfy the various constraints. No acceptable cases.
    
   OUT=nan; 
end


end
