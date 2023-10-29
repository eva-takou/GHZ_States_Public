%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 29, 2023
%--------------------------------------------------------------------------

%%
path = '/Users/evatakou/Documents/MATLAB/GHZ_States_Public/Simulations/GHZ_Data_of_Sequential_FINAL/GHZ3_Sequential.mat';
load(path,"OUT")
clc

[s0,s1,wL]  = load_fixed_params;
s00         = [1 0 ; 0 0];
s11         = [0 0 ; 0 1];
ket0        = [1;0];   


CaseNum =  3;

%%

At = OUT.A_Target{CaseNum};
Bt = OUT.B_Target{CaseNum};
Au = OUT.A_Unwanted{CaseNum};
Bu = OUT.B_Unwanted{CaseNum};

Au = Au(1:4); %Just consider a few of the unwanted nuclei.
Bu = Bu(1:4);

times = OUT.Opt_Unit_Times{CaseNum};
iters = OUT.Opt_Iters{CaseNum};

Ntarget     = length(At);

[rho_red_Analy,~,~,~,PSI0]=reduced_density_matrix_CR_type(At,Bt,Au,Bu,times,iters);

%Test if we found the correct reduced density matrix
k    = 1;
test = SuperClass_Sequences(wL,[At,Au],[Bt,Bu],s0,s1,length(At)+length(Au),k,1);
test = test.CPMG(times,iters);
U    = test.Uevol;

clc

for jj=1:length(Au)
   
    PSI0 = kron(PSI0,ket0);
end

PSI  = U*PSI0;
rho  = PSI*PSI';
NQ   = log2(length(PSI));
Nunw = length(Au);

rho_red = ptrace(rho,...
                 Ntarget+2:Ntarget+1+Nunw,...
                 repmat(2,[1,NQ]));

              
test_cond = norm(rho_red-rho_red_Analy);

if abs(test_cond)>1e-9
    
    error('Wrong evaluation of reduced density matrix.')
    
else
    
    disp('Evaluation of reduced density matrix is correct.')
    
    
end




