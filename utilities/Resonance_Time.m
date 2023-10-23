function tres = Resonance_Time(k,Nuclear_Indx,wL,A,B,s0,s1)
%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
%Resonance time.
%Input: k: index of resonance
%       INDX: indx of particular nuclear spin (1-27 from Delft)
%       wL: Larmor frequency
%       A, B: Parallel, perpendicular HF components of nuclei
%       s0,s1: Electron's spin projections

c = 2*pi*1e-3; %Conversion factor to \omega (rad*MHz)

W0 = @(A,B)  sqrt( (wL+s0*A)^2 +(s0*B)^2)*c;
W1 = @(A,B)  sqrt( (wL+s1*A)^2 +(s1*B)^2)*c;
W  = @(A,B)  W0(A,B)+W1(A,B);

time = @(A,B) 4*pi*(2*k-1)/(W(A,B));

tres = time(A(Nuclear_Indx),B(Nuclear_Indx));


end
