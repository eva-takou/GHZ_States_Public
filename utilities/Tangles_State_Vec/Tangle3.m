function tau=Tangle3(psi)
%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 29, 2023
%--------------------------------------------------------------------------
%Calculate the 3-tangle of a three-qubit state. 
%This can be done using the method of Ref:
%V. Coffman et al, Phys. Rev. A 61, 052306 (2000).

if log2(length(psi))~=3
    
    error('The input state is not a three-qubit state.')
    
end


psi = psi/norm(psi);

c000 = psi(1); c001 = psi(2); c010 = psi(3); c011 = psi(4); 
c100 = psi(5); c101 = psi(6); c110 = psi(7); c111 = psi(8);


d1 = c000^2*c111^2+c001^2*c110^2 + c010^2*c101^2 + c100^2*c011^2;

d2 = c000*c111*c011*c100 +...
     c000*c111*c101*c010 +...
     c000*c111*c110*c001 +...
     c011*c100*c101*c010 +...
     c011*c100*c110*c001 +...
     c101*c010*c110*c001;
 
d3 = c000 * c110 * c101 * c011 + c111 * c001 * c010 * c100;


tau=4*abs(d1-2*d2+4*d3);




 

end