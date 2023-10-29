%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
% Script to test numerical with analytical evaluation of non-unitary
% M-tangling power
clc;
clear
close all;

%% Set basic parameters

Compositions = 2;
wL           = 314;
times        = 1+(10-1)*rand(1,Compositions);    
iters        = randi(100,[1,Compositions]);
s0           = 0;      
s1           = -1;

Nnuc = 6;

a = 20;
b = 200;

A =  a + (b-a).*rand(Nnuc,1);
B =  a + (b-a).*rand(Nnuc,1);



%% Target nuclei

LT        = 3;
Lunw      = Nnuc-LT;
M         = LT+1;
d         = 2;
indx_Targ = 1:LT;
indx_unw  = setxor(1:Nnuc,indx_Targ);

Aunw = A(indx_unw);
Bunw = B(indx_unw);
At   = A(indx_Targ);
Bt   = B(indx_Targ);

%% Get the evolution of target/unwanted spins

for jj=1:Lunw
   
    spin=SubClass_U4Operations(wL,Aunw(jj),Bunw(jj),s0,s1,1,1,1);
    spin=spin.CPMG(times,iters);
    spin=spin.Rot_Angles_And_Axes;
    phi0_u(jj) = spin.angles{1};
    phi1_u(jj) = spin.angles{2};
    n0_u(jj,:) = spin.axes{1};
    n1_u(jj,:) = spin.axes{2};
end


for jj=1:LT
   
    spin=SubClass_U4Operations(wL,At(jj),Bt(jj),s0,s1,1,1,1);
    spin=spin.CPMG(times,iters);
    spin=spin.Rot_Angles_And_Axes;
    phi0_t(jj) = spin.angles{1};
    phi1_t(jj) = spin.angles{2};
    n0_t(jj,:) = spin.axes{1};
    n1_t(jj,:) = spin.axes{2};
end
%% Get evolution operator of target space

temp=SuperClass_Sequences(wL,At,Bt,s0,s1,length(At),1,1);
temp=temp.CPMG(times,iters);
U0 =temp.Uevol;
%% Get the evolution of total space

temp=SuperClass_Sequences(wL,A,B,s0,s1,length(A),1,1);
temp=temp.CPMG(times,iters);
U =temp.Uevol;

%% Get the Kraus operators
temp=SubClass_Ent_and_Fid;
Ek=temp.Get_Kraus(U,Nnuc,indx_Targ);
%% Evaluate analytical expression of M-tangling power
temp=SubClass_Ent_and_Fid;
temp=temp.NonUniMwayEp_CR_Analytical(Nnuc,length(indx_Targ),phi0_u,phi1_u,n0_u,n1_u,phi0_t,phi1_t,n0_t,n1_t);

non_uni_epM = temp.epM_nonuni;


%% Now evaluate in purely numerical way

temp=0;
Pm=1;
Pp=1;

if (-1)^M==1 %even
   
    for ii=1:M
        Pm=Pm*1/2*(eye(2^(2*M))-ArbDimSWAP(ii,ii+M,2*M));
    end
 
    for ii=1:M
        Pp=Pp*1/2*(eye(2^(2*M))+ArbDimSWAP(ii,ii+M,2*M));
    end
    
    
else %odd
    
    for ii=2:M
        Pm=Pm*1/2*(eye(2^(2*M))-ArbDimSWAP(ii,ii+M,2*M));
    end

    for ii=1:(LT+1)
        Pp=Pp*1/2*(eye(2^(2*M))+ArbDimSWAP(ii,ii+M,2*M));
    end
    
end


Omp0 = 1/(d+1)^M*Pp;

for r=1:size(Ek,3)
    
    for s=1:size(Ek,3)
    
        Er   = Ek(:,:,r);
        Es   = Ek(:,:,s);
        
        temp = temp + trace(kron(Er,Es)*Omp0*kron(Er',Es')*Pm);
    
    end
    
end

temp = 2^M*temp;

%%

abs(temp-non_uni_epM)

