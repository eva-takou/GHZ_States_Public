clear
M = 4;  
d = 2;  
path_Sequential='/Users/evatakou/Documents/MATLAB/Defect_Nuclear_GHZ_states/GHZ_states_Preparation/Simulations/GHZ_Data_of_Sequential_NEW/GHZ4_Sequential.mat';
load(path_Sequential,"OUT");
S = OUT;


%% 1e-3

tol=5e-3;

ii=1;

Target_Nuclei   = S.Target_Nuclei{ii};
Unwanted_Nuclei = S.Unwanted_Nuclei_Names{ii};
temp            = S.EP_Unwanted{ii};
indices         = temp>tol; 
Unwanted_Nuclei = Unwanted_Nuclei(indices);
Ntarget         = length(Target_Nuclei);
Nnuc            = length(Target_Nuclei)+length(Unwanted_Nuclei);

%Parameters of unwanted spins evolution
phi0 = S.phi0N{ii}(indices);
phi1 = S.phi1N{ii}(indices);
n0   = S.n0{ii}(indices,:);
n1   = S.n1{ii}(indices,:);

%Parameters of target spins evolution
th0  = S.phi0N{ii}(Target_Nuclei);
th1  = S.phi1N{ii}(Target_Nuclei);
u0   = S.n0{ii}(Target_Nuclei,:);
u1   = S.n1{ii}(Target_Nuclei,:);

temp     = SubClass_Ent_and_Fid;
temp     = temp.NonUniMwayEp_CR(Nnuc,Ntarget,phi0,phi1,n0,n1,th0,th1,u0,u1);
ep_E(ii) = temp.epM_nonuni/(d/(d+1))^M; %Scale by max value
ep_U(ii) = prod(S.EP_Target{ii}); 


clc
ep_E;
ep_U;


G1 =@(n0,n1,phi0,phi1) (cos(phi0/2)*cos(phi1/2)+dot(n0,n1)*sin(phi0/2)*sin(phi1/2))^2;
g1_prod = 1;

for jj=1:length(Unwanted_Nuclei)
    
    g1_prod = g1_prod * G1(n0(jj,:),n1(jj,:),phi0(jj),phi1(jj));
    
end

guess(ii)=(1+g1_prod)/2*prod(S.EP_Target{ii});

disp(['Guess expression:',num2str(guess)])
disp(['Numerical expression:',num2str(ep_E)])


%% 5e-4
clc
ii=3;

Target_Nuclei   = S.Target_Nuclei{ii};
Unwanted_Nuclei = S.Unwanted_Nuclei_Names{ii};
temp            = S.EP_Unwanted{ii};
indices         = temp>1e-2; 
Unwanted_Nuclei = Unwanted_Nuclei(indices);
Ntarget         = length(Target_Nuclei);
Nnuc            = length(Target_Nuclei)+length(Unwanted_Nuclei);

%Parameters of unwanted spins evolution
phi0 = S.phi0N{ii}(Unwanted_Nuclei);
phi1 = S.phi1N{ii}(Unwanted_Nuclei);
n0   = S.n0{ii}(Unwanted_Nuclei,:);
n1   = S.n1{ii}(Unwanted_Nuclei,:);

%Parameters of target spins evolution
th0  = S.phi0N{ii}(Target_Nuclei);
th1  = S.phi1N{ii}(Target_Nuclei);
u0   = S.n0{ii}(Target_Nuclei,:);
u1   = S.n1{ii}(Target_Nuclei,:);

temp     = SubClass_Ent_and_Fid;
temp     = temp.NonUniMwayEp_CR(Nnuc,Ntarget,phi0,phi1,n0,n1,th0,th1,u0,u1);
temp2     = temp.NonUniMwayEp_CR_Faster(Nnuc,Ntarget,phi0,phi1,n0,n1,th0,th1,u0,u1);
ep_E(ii) = temp.epM_nonuni/(d/(d+1))^M; %Scale by max value
ep_E2(ii) = temp2.epM_nonuni/(d/(d+1))^M; %Scale by max value
ep_U(ii) = prod(S.EP_Target{ii}); 



ep_E %0.8463
ep_U
ep_E2

%% 4e-4

ii=1;

Target_Nuclei   = S.Target_Nuclei{ii};
Unwanted_Nuclei = S.Unwanted_Nuclei_Names{ii};
temp            = S.EP_Unwanted{ii};
indices         = temp>4e-4; 
Unwanted_Nuclei = Unwanted_Nuclei(indices);
Ntarget         = length(Target_Nuclei);
Nnuc            = length(Target_Nuclei)+length(Unwanted_Nuclei);

%Parameters of unwanted spins evolution
phi0 = S.phi0N{ii}(Unwanted_Nuclei);
phi1 = S.phi1N{ii}(Unwanted_Nuclei);
n0   = S.n0{ii}(Unwanted_Nuclei,:);
n1   = S.n1{ii}(Unwanted_Nuclei,:);

%Parameters of target spins evolution
th0  = S.phi0N{ii}(Target_Nuclei);
th1  = S.phi1N{ii}(Target_Nuclei);
u0   = S.n0{ii}(Target_Nuclei,:);
u1   = S.n1{ii}(Target_Nuclei,:);

temp     = SubClass_Ent_and_Fid;
temp     = temp.NonUniMwayEp_CR(Nnuc,Ntarget,phi0,phi1,n0,n1,th0,th1,u0,u1);
ep_E(ii) = temp.epM_nonuni/(d/(d+1))^M; %Scale by max value
ep_U(ii) = prod(S.EP_Target{ii}); 

ep_E %0.8463
ep_U


%% 1e-2

ii=1;

Target_Nuclei   = S.Target_Nuclei{ii};
Unwanted_Nuclei = S.Unwanted_Nuclei_Names{ii};
temp            = S.EP_Unwanted{ii};
indices         = temp>1e-2; 
Unwanted_Nuclei = Unwanted_Nuclei(indices);
Ntarget         = length(Target_Nuclei);
Nnuc            = length(Target_Nuclei)+length(Unwanted_Nuclei);

%Parameters of unwanted spins evolution
phi0 = S.phi0N{ii}(Unwanted_Nuclei);
phi1 = S.phi1N{ii}(Unwanted_Nuclei);
n0   = S.n0{ii}(Unwanted_Nuclei,:);
n1   = S.n1{ii}(Unwanted_Nuclei,:);

%Parameters of target spins evolution
th0  = S.phi0N{ii}(Target_Nuclei);
th1  = S.phi1N{ii}(Target_Nuclei);
u0   = S.n0{ii}(Target_Nuclei,:);
u1   = S.n1{ii}(Target_Nuclei,:);

temp     = SubClass_Ent_and_Fid;
temp     = temp.NonUniMwayEp_CR(Nnuc,Ntarget,phi0,phi1,n0,n1,th0,th1,u0,u1);
ep_E(ii) = temp.epM_nonuni/(d/(d+1))^M; %Scale by max value
ep_U(ii) = prod(S.EP_Target{ii}); 

ep_E %0.8503
ep_U
