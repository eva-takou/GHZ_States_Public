function OUT=minimize_MTangle(At,Bt,Aunw,Bunw,times,iters)
%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
%
%Function to minimize the M-tangle of a mixed state of rank 2 via convex hull
%construction.[See Ref. https://dx.doi.org/10.1103/PhysRevA.77.032310]
%
%Input: At: Array of parallel HF components of target nuclei
%       Bt: Array of perpendicular HF components of target nuclei
%       Aunw: Array of parallel HF components of unwanted nuclei 
%       Bunw: Array of perpendicular HF components of unwanted nuclei 
%       times: unit times of composite sequences of sequential scheme
%       iters: iterations of each unit of the composite sequences of sequential scheme
%
%Output: OUT: A structwith various fields that contains information
%             about the mixed state
%
%--------------------------------------------------------------------------

[~,evecs,evals,PSI_Pure]=reduced_density_matrix_CR_type(At,Bt,Aunw,Bunw,times,iters);

evecs = evecs(:,1:2);
evals = evals(1:2);

psi_trial = @(gamma) sqrt(evals(1))*evecs(:,1)-exp(1i*gamma)*sqrt(1-evals(1))*evecs(:,2);

NumQubits   = log2(length(evecs(:,1)));
gamma_range = 0:0.005*pi:2*pi;
tau         = zeros(1,length(gamma_range));

for ii = 1:length(gamma_range)
    
    if NumQubits ==3
        
        tau(ii) =Tangle3(psi_trial(gamma_range(ii)));
        
    elseif (-1)^NumQubits ==1
        
        tau(ii) =TangleEven(psi_trial(gamma_range(ii)));
        
    elseif NumQubits == 5
        
        tau(ii) =Tangle5(psi_trial(gamma_range(ii)));
        
    end
    
end

[minVal,indx]=min(tau);

OUT.gamma_Opt = gamma_range(indx);
OUT.minTau    = minVal;
OUT.p         = evals(1);
OUT.tau_Vec1  = Tangle3(evecs(:,1));
OUT.tau_Vec2  = Tangle3(evecs(:,2));

%Find the tangle of the Pure_State
    
if NumQubits==3
    
    OUT.PureTangle = Tangle3(PSI_Pure);

elseif (-1)^NumQubits==1

    OUT.PureTangle = TangleEven(PSI_Pure);

elseif NumQubits==5
    
    OUT.PureTangle = Tangle5(PSI_Pure);

end
    



end





