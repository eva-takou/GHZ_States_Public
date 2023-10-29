function tau=TangleEven(psi)
%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 29, 2023
%--------------------------------------------------------------------------
%Calculate the M-tangle of any state with even # of qubits.
%This can be done using the method of Ref:
%V. Coffman et al, Phys. Rev. A 61, 052306 (2000).
%--------------------------------------------------------------------------

L = length(psi);
L = log2(L); 

if (-1)^L~=1
   
    error('The input state does not correspond to an even # of qubits.')
    
end

psi = psi/norm(psi);     

Y = [0 -1i ; 1i 0];    
Yall = 1;

for ii=1:L
    
    Yall = kron(Yall,Y);
    
end

tau = abs(  psi' * Yall * conj(psi)   )^2;



end
