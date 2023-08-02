function tau=TangleEven(psi)

L = length(psi);
L = log2(L); 

psi = psi/norm(psi);     

Y = [0 -1i ; 1i 0];    
Yall = 1;

for ii=1:L
    
    Yall = kron(Yall,Y);
    
end

tau = abs(  psi' * Yall * conj(psi)   )^2;



end
