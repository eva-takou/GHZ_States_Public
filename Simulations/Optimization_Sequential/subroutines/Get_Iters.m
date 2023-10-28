function [t_opt,N]=Get_Iters(Spin_Indx,Res,wL,A,B,s0,s1,Number_Of_Maxima,Time_Max,Unwanted_One_Tangle_Tol)
%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 27, 2023
%--------------------------------------------------------------------------
%Output all possible Iters and 1 optimal time of a particular resonance and
%nuclear spin such that the unwanted one-tangles are below the tolerance.
%Input: Spin_Indx \in [1,27]
%       Res: resonance index
%       wL,A,B: Larmor frequency and HF parameters of 27 nuclei
%       s0,s1: Electron's projections
%       Number_Of_Maxima: how many maxima of target one-tangle to keep 
%       Time_Max: Time constraint of sequence
%       Unwanted_One_Tangle_Tol: Tolerance of unwanted one-tangles.
%Output: optimal unit time and iterations N.


Nnuc    = 27;
TIME    = Resonance_Time(Res,Spin_Indx,wL,A,B,s0,s1);
tt      = TIME-0.2:5e-3:TIME+0.2; %Vary the time by \pm 0.2 \mu s
numSpin = 1;
n0n1    = zeros(1,length(tt));

for ii=1:length(tt) %Find when the dot product==-1
    
NucSpin = SubClass_U4Operations(wL,A(Spin_Indx),B(Spin_Indx),s0,s1,numSpin,Res,1);  %1 iteration
NucSpin  = NucSpin.CPMG(tt(ii));
NucSpin  = NucSpin.Rot_Angles_And_Axes;
n0n1(ii) = dot(NucSpin.axes{1},NucSpin.axes{2} );

end

[~,indx]=min(n0n1);    t_opt = tt(indx);
clear tt
clear n0n1
NucSpin = SubClass_U4Operations(wL,A(Spin_Indx),B(Spin_Indx),s0,s1,numSpin,Res,1);  %1 iteration)
NucSpin = NucSpin.CPMG(t_opt);
NucSpin = NucSpin.Expected_Maxima(Number_Of_Maxima);
N       = NucSpin.Nmax;

for ii=1:length(N)
    
   if ~isreal(N(ii)) || N(ii)<0 
   
       error('Error in the evaluation of iterations that maximize one-tangles.')
       
   end
    
end

N = N(N*t_opt<Time_Max); %Restrict N such that we are within time limits.

if isempty(N) 
    
    t_opt = nan;
    N     = nan;
    
    return
    
end

%Remove N where epunwanted>unwanted one-tangle tolerance 
for ii=1:length(N)
    
    spin=1;
    
    while spin<=Nnuc
         
        NucSpin = SubClass_U4Operations(wL,A(spin),B(spin),s0,s1,numSpin,Res,N(ii));  
        NucSpin = NucSpin.CPMG(t_opt);
        NucSpin = NucSpin.Makhlin_Inv;
        Ep      = NucSpin.Ep/(2/9);
        
        if spin~=Spin_Indx && Ep>Unwanted_One_Tangle_Tol
            
            N(ii)=nan;
            break
            
        end
        spin=spin+1;
    end
    
end

N = N(~isnan(N));

if isempty(N)
    
    t_opt = nan;
    N     = nan;
end

                   
end
