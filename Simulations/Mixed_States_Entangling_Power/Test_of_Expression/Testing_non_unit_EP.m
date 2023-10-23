%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------


%% First, test for a single unwanted spin, 2 target ones

clearvars; clc;

s00  =[1 0 ; 0 0 ]; s11  =[0 0 ; 0 1 ];
X    =[0 1; 1 0]; 
Y    =[0, -1i ; 1i 0]; 
Z    =[1 0 ; 0 -1];

Rn = @(phi,n) expm(-1i*phi/2*(n(1)*X+n(2)*Y+n(3)*Z));

%% Set basic parameters

Compositions = 2;
wL           = 314;
times        = 1+(10-1)*rand(1,Compositions);    
iters        = randi(100,[1,Compositions]);
s0           = 0;      
s1           = -1;

%%

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

%%
N=1;
option = "Only_Kraus";
%total system
Test         = SuperClass_Sequences(wL,A,B,s0,s1,Nnuc,1,N);
Test         = Test.CPMG(times,iters);
U            = Test.Uevol;

%just target subspace
Test         = SuperClass_Sequences(wL,A(1:LT),B(1:LT),s0,s1,LT,1,N);
Test         = Test.CPMG(times,iters);
U0           = Test.Uevol;

%just electron + single spin subspace

Test         = SubClass_U4Operations(wL,A(1),B(1),s0,s1,1,1,N);
Test         = Test.CPMG(times,iters);
U0_en        = Test.Uevol;

%% Get the Kraus for the system of electron + single spin
%  or electron+ multiple target spins

Empty        = SubClass_Ent_and_Fid;
Empty        = Empty.GKFnew(U,{U0_en,U0},Nnuc,indx_Targ,option);

[~,~,LKraus] = size(Empty.Kraus_MultiQ);

for ii=1:Nnuc
    Test         = SubClass_U4Operations(wL,A(ii),B(ii),s0,s1,1,1,N);
    Test         = Test.CPMG(times,iters);
    Test         = Test.Rot_Angles_And_Axes;
    Test         = Test.Makhlin_Inv;
    ep(ii)       = Test.Ep/(2/9);
    phi0(ii)     = Test.angles{1};
    phi1(ii)     = Test.angles{2};
    n0(ii,:)     = Test.axes{1};
    n1(ii,:)     = Test.axes{2};
end

%% Calculate the entangling power based on full numerical expression

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

for r=1:LKraus
    
    for s=1:LKraus
    
        Er   = Empty.Kraus_MultiQ(:,:,r);
        Es   = Empty.Kraus_MultiQ(:,:,s);
        
        temp = temp + trace(kron(Er,Es)*Omp0*kron(Er',Es')*Pm);
    
    end
    
end

temp = 2^M*temp;

disp(['Full numerical M-way EP: ',num2str(real(temp))])

%% Now use the expression we have in the SubClass_Ent_and_Fid

I1  = indx_unw;
% 
% 
% c0_el0 =@(indx) cos(phi0(indx)/2)-1i*n0(indx,3)*sin(phi0(indx)/2);
% c1_el0 =@(indx) -1i*sin(phi0(indx)/2)*(n0(indx,1)+1i*n0(indx,2));
% 
% c0_el1 =@(indx) cos(phi1(indx)/2)-1i*n1(indx,3)*sin(phi1(indx)/2);
% c1_el1 =@(indx) -1i*sin(phi1(indx)/2)*(n1(indx,1)+1i*n1(indx,2));
% 
% G1   =@(indx)  ( cos(phi0(indx)/2)^2 + dot(n0(indx,:),n1(indx,:))*sin(phi0(indx)/2)^2 )^2;

%Kraus 1: nucleus in 0 electron in 0/1
% 
% f0_Kraus1 = c0_el0(I1);
% f1_Kraus1 = c0_el1(I1);
% 
% %Kraus 2: nucleus in 1 electron in 0/1
% 
% f0_Kraus2 = c1_el0(I1);
% f1_Kraus2 = c1_el1(I1);


% prefactor =  (d/(d+1))^M * (1-G1(1))*(1-G1(2));
% 
% 
% test2 = prefactor*1/2*( abs(f0_Kraus1)^2*abs(f1_Kraus2)^2 + abs(f0_Kraus2)^2*abs(f1_Kraus1)^2 +...
%                         abs(f0_Kraus2)^2*abs(f1_Kraus1)^2 + abs(f0_Kraus1)^2*abs(f1_Kraus2)^2 +...
%                         abs(f0_Kraus1)^2*abs(f1_Kraus1)^2 + abs(f0_Kraus1)^2*abs(f1_Kraus1)^2 +...
%                         abs(f0_Kraus2)^2*abs(f1_Kraus2)^2 + abs(f0_Kraus2)^2*abs(f1_Kraus2)^2);
% 

% disp(['Numerical M-way EP: ',num2str(test2)])
% 
% %let's reconstruct the Kraus:
% 
% e0 = [1;0]; e0 = kron(eye(2^(LT+1)),e0);
% e1 = [0;1]; e1 = kron(eye(2^(LT+1)),e1);
% 
% 
% Ek1 = e0'*U*e0; Ek2 = e1'*U*e0; 
% 
% Ek1_Alt = f0_Kraus1*kron(s00,kron(Rn(phi0(1),n0(1,:)),Rn(phi0(2),n0(2,:))))+...
%           f1_Kraus1*kron(s11,kron(Rn(phi1(1),n1(1,:)),Rn(phi1(2),n1(2,:))));
% 
% Ek2_Alt = f0_Kraus2*kron(s00,kron(Rn(phi0(1),n0(1,:)),Rn(phi0(2),n0(2,:))))+...
%           f1_Kraus2*kron(s11,kron(Rn(phi1(1),n1(1,:)),Rn(phi1(2),n1(2,:))));
%        

% disp(['deviation of Kraus 1 from numerical:',num2str(norm(Empty.Kraus_MultiQ(:,:,1)-Ek1))])
% disp(['deviation of Kraus 2 from numerical:',num2str(norm(Empty.Kraus_MultiQ(:,:,2)-Ek2))])
% 
% %up to global \pm sign
% disp('============================')
% disp(['deviation of Alt Kraus 1 from numerical:',num2str(norm(Empty.Kraus_MultiQ(:,:,1)+Ek1_Alt))])
% disp(['deviation of Alt Kraus 1 from numerical:',num2str(norm(Empty.Kraus_MultiQ(:,:,2)+Ek2_Alt))])
% disp(['deviation of Alt Kraus 1 from numerical:',num2str(norm(Empty.Kraus_MultiQ(:,:,1)-Ek1_Alt))])
% disp(['deviation of Alt Kraus 1 from numerical:',num2str(norm(Empty.Kraus_MultiQ(:,:,2)-Ek2_Alt))])

Empty=SubClass_Ent_and_Fid;

test2=Empty.NonUniMwayEp_CR(Nnuc,LT,phi0(I1),phi1(I1),n0(I1,:),n1(I1,:),phi0(1:LT),phi1(1:LT),...
                     n0(1:LT,:),n1(1:LT,:));

disp(['Evaluation from Class:',num2str(test2.epM_nonuni)])

%% Now use the guess expression

G1     = @(phi0,phi1,n0,n1) (cos(phi0/2)*cos(phi1/2)+dot(n0,n1)*sin(phi0/2)*sin(phi1/2))^2;
g1prod = 1;

for jj=1:Lunw
    
   g1prod = g1prod*(G1(phi0(I1(jj)),phi1(I1(jj)),n0(I1(jj),:),n1(I1(jj),:)));
    
end

g1T=1;

for jj=1:LT
    
    g1T = g1T*(1-G1(phi0(jj),phi1(jj),n0(jj,:),n1(jj,:)));
end


guess=(d/(d+1))^M*(1+g1prod)/2*g1T;

disp(['Guess expr=',num2str(guess)])


