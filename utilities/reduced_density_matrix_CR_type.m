function [rho_red_Analy,evecs,evals,PSI_Pure,PSI_0]=reduced_density_matrix_CR_type(At,Bt,Aunw,Bunw,times,iters)
%Get the reduced density matrix after tracing out unwanted nuclei.
%Assume the CPMG sequence. For composite, we are composing CPMG sequences
%of different unit times and iterations.
%
%Input: Target HF parameters (At and Bt)
%       Unwanted HF parameters (Au and Bu)
%       Unit time(s) of CPMG sequences  (composite or just 1 CPMG sequence)
%       Iterations(s) of CPMG sequences (composite or just 1 CPMG sequence)
%Output: rho_red_Analy: reduced density matrix
%        evecs: eigenvectors of the mixed state
%        evals: eigenvalues corresponding to the eigenvectors
%        PSI_Pure: pure state, assuming the evolution of only the target
%        subspace (consisting of the electron and target nuclei)
%        rho_el0: Initial state of target subspace
if length(At)~=2 || length(Bt)~=2
   
    error('This script is for 2 nuclear spins and the electron only. Consider generalizing this.')
    
end

%%========= Params =======================================================

[s0,s1,wL]  = load_fixed_params;
s00         = [1 0 ; 0 0];
s11         = [0 0 ; 0 1];
ket0        = [1;0];   
rho         = @(psi) psi*psi';
Ntarget     = length(At);

%%=========  Get the evolution of all target spins ========================

R0_T = cell(1,Ntarget);
R1_T = cell(1,Ntarget);

for ii=1:Ntarget
    
    Spin     = SubClass_U4Operations(wL,At(ii),Bt(ii),s0,s1,1,1,1);
    Spin     = Spin.CPMG(times,iters);
    R0_T{ii} = Spin.Uevol(1:2,1:2);  
    R1_T{ii} = Spin.Uevol(3:4,3:4);
    
end

r0 = 1; r1=1;

for ii=1:Ntarget
   
    r0 = kron(r0,R0_T{ii});
    r1 = kron(r1,R1_T{ii});
    
end

CR                     = kron(s00,r0)+kron(s11,r1);
[psi_el,psi_n1,psi_n2] = Construct_Init_State(CR); %Initial states of target space
psi_nuc_Unw            = ket0;                     %Assume unwanted spins in |0> state
rho_ele0               = rho(psi_el);

%%======== Evolve unwanted nuclei =========================================

F00 = 1; F11 = 1; F01 = 1; F10 = 1; %f00 = 1 = f11 always

Nunw = length(Aunw);

for ii=1:Nunw 
    
    Spin = SubClass_U4Operations(wL,Aunw(ii),Bunw(ii),s0,s1,1,1,1);
    Spin = Spin.CPMG(times,iters); %Here we can change the sequence
    
    Rn0  = Spin.Uevol(1:2,1:2);
    Rn1  = Spin.Uevol(3:4,3:4);
    
    F00 = F00 * trace(Rn0*rho(psi_nuc_Unw)*Rn0'); %=1
    F11 = F11 * trace(Rn1*rho(psi_nuc_Unw)*Rn1'); %=1
    F01 = F01 * trace(Rn0*rho(psi_nuc_Unw)*Rn1');
    F10 = F10 * trace(Rn1*rho(psi_nuc_Unw)*Rn0'); 
    
end

%%======= Evolve the target nuclei ========================================

rho_nuc_T_00 = R0_T{1}*rho(psi_n1)*R0_T{1}';
rho_nuc_T_01 = R0_T{1}*rho(psi_n1)*R1_T{1}';
rho_nuc_T_10 = R1_T{1}*rho(psi_n1)*R0_T{1}';
rho_nuc_T_11 = R1_T{1}*rho(psi_n1)*R1_T{1}';

for ii=2:Ntarget
    
    rho_nuc_T_00 = kron(rho_nuc_T_00,R0_T{ii}*rho(psi_n2)* R0_T{ii}');
    rho_nuc_T_01 = kron(rho_nuc_T_01,R0_T{ii}*rho(psi_n2)* R1_T{ii}');
    rho_nuc_T_10 = kron(rho_nuc_T_10,R1_T{ii}*rho(psi_n2)* R0_T{ii}');
    rho_nuc_T_11 = kron(rho_nuc_T_11,R1_T{ii}*rho(psi_n2)* R1_T{ii}');
    
end

PSI_0    = kron(psi_el,kron(psi_n1,psi_n2));
PSI_Pure = CR*PSI_0;


%======== Construct the analytical reduced density matrix ================

rho_red_Analy = F00*kron(s00*rho_ele0*s00,rho_nuc_T_00) +...
                F01*kron(s00*rho_ele0*s11,rho_nuc_T_01) +...
                F10*kron(s11*rho_ele0*s00,rho_nuc_T_10) +...
                F11*kron(s11*rho_ele0*s11,rho_nuc_T_11);

[evecs,evals]=sorted_evals_evecs(rho_red_Analy); %get evals and evecs

if norm(evecs(:,1)'*evecs(:,2))>1e-6
    
   error('The e-evecs are not orthogonal.') 
   
end

if evals(2)<1e-8
    
    error('Second eval=0.')
    
end

disp(['Tangle of 1st evector:',num2str(Tangle3(evecs(:,1)))])

disp(['Tangle of 2nd evector:',num2str(Tangle3(evecs(:,2)))])
disp('=========================================')

end

function [psi1,psi2,psi3]=Construct_Init_State(CR)
%To optimize the initial state only for 3 qubits.

ket0 = [1;0];
ket1 = [0;1];
ketP = 1/sqrt(2)*(ket0+ket1);

TOL      = 0.95; 
psi_test = kron(ketP,kron(ket0,ket0));

if Tangle3(CR*psi_test)>TOL
    
    psi1 = ketP;
    psi2 = ket0;
    psi3 = ket0;
    
    return
else
    
    disp('ENTERED OPTIMIZATION OF INITIAL STATE')
    
end

psiQ = @(th,ph) cos(th/2)*ket0 + exp(1i*ph)*sin(th/2)*ket1;
step = 0.1*pi;
th   = 0:step/2:pi;
ph   = 0:step:2*pi;
cnt  = 0;

for i1=1:length(th)
    
    for j1=1:length(ph)
        
        for i2=1:length(th)
            for j2=1:length(ph)
                
                for i3=1:length(th)
                    for j3=1:length(ph)
                        
                     psi1 = psiQ(th(i1),ph(j1));
                     psi2 = psiQ(th(i2),ph(j2));
                     psi3 = psiQ(th(i3),ph(j3));
                     psi  = kron(psi1,kron(psi2,psi3));
                     cnt  = cnt+1;
                     
                     tau3(i1,j1,i2,j2,i3,j3) = Tangle3(CR*psi);
                        
                    end
                end
            end
        end
    end
    
end

[~,indx]            = max(tau3,[],'all','linear');
[I1,J1,I2,J2,I3,J3] = ind2sub(size(tau3),indx);

psi1 = psiQ(th(I1),ph(J1));
psi2 = psiQ(th(I2),ph(J2));
psi3 = psiQ(th(I3),ph(J3));

end
