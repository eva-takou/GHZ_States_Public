%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 29, 2023
%--------------------------------------------------------------------------
clear;
clc
%% Test the 3-tangle of GHZ

ket0=[1;0];
ket1=[0;1];

GHZsize=3;

psiA=ket0;
psiB=ket1;

for jj=1:GHZsize-1
    
    psiA = kron(psiA,ket0);
    psiB = kron(psiB,ket1);
    
end

Psi = psiA+psiB;
Psi = Psi/norm(Psi);

Tangle3(Psi)
%% Test the 5-tangle of GHZ

GHZsize=5;

psiA=ket0;
psiB=ket1;

for jj=1:GHZsize-1
    
    psiA = kron(psiA,ket0);
    psiB = kron(psiB,ket1);
    
end

Psi = psiA+psiB;
Psi = Psi/norm(Psi);

Tangle5(Psi)

%% Test the 3-Tangle of state equivalent to GHZ
X=[0 1;1 0];
Y=[0 -1i;1i 0];
Z=[1 0 ; 0 -1];

Rn=@(a,n) expm(-1i*a/2*(n(1)*X+n(2)*Y+n(3)*Z));


GHZsize=3;

psiA=1;
psiB=1;

for jj=1:GHZsize
       
    psiA=kron(psiA,ket0);
    psiB=kron(psiB,ket1);
    
    
end

Psi = psiA+psiB;

ph=2*pi*rand(1);
th=pi*rand(1);
a=2*pi*rand(1);


nx = sin(th)*cos(ph);
ny = sin(th)*sin(ph);
nz = cos(th);

Rot = Rn(a,[nx,ny,nz]);



Psi = kron(Rot,eye(4))*Psi;
Psi = Psi/norm(Psi);

Tangle3(Psi)

%% Test the 5 tangle of GHZ-equivalent state

GHZsize=5;

psiA=1;
psiB=1;

for jj=1:GHZsize
       
    psiA=kron(psiA,ket0);
    psiB=kron(psiB,ket1);
    
end

Psi = psiA+psiB;

ph=2*pi*rand(1);
th=pi*rand(1);
a=2*pi*rand(1);


nx = sin(th)*cos(ph);
ny = sin(th)*sin(ph);
nz = cos(th);

Rot = Rn(a,[nx,ny,nz]);

Psi = kron(Rot,eye(2^4))*Psi;
Psi = Psi/norm(Psi);

Tangle5(Psi)
%% Test the 3-tangle of W state


W=kron(ket0,kron(ket0,ket1))+kron(ket1,kron(ket0,ket0))+kron(ket0,kron(ket1,ket0));

W=W/norm(W);

Tangle3(W)










