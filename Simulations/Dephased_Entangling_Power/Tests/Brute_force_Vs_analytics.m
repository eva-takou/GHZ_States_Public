%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------


clear

clc

%%

s0=0;
s1=-1;
c     = 2*pi*1e-3;
A1=60;
B1=30;
A2=70;
B2=40;
A=[A1,A2];
B=[B1,B2];

wl=314;

Om0 = sqrt((wl+s0*A1)^2+(s0*B1)^2);
Om1 = sqrt((wl+s1*A1)^2+(s1*B1)^2);

k     = 1;
t     = 4*pi*(2*k-1)/(Om0+Om1)*1/c;
t     = [t,0.9*t];

N     = [1,1];
theta = 1/400;
l0    = 0.4;
l1    = 0.6;

%%

for jj=1:length(A)

    temp=SubClass_U4Operations(wl,A(jj),B(jj),s0,s1,1,1,1);
    temp=temp.CPMG(t,N);
    temp=temp.Rot_Angles_And_Axes;
    phi0(jj)=temp.angles{1};
    phi1(jj)=temp.angles{2};
    n0{jj}=temp.axes{1};
    n1{jj}=temp.axes{2};
end

%%
M=length(A)+1;
d=2;
clc

temp     = SubClass_Ent_and_Fid;
temp     = temp.dephased_epM_CR_Corr(n0,n1,phi0,phi1,l0,l1,theta,t,N);
epM      = temp.epM_dephased/(d/(d+1))^M

%%

temp   = SubClass_Ent_and_Fid;
temp   = temp.dephased_epM_brute_force(wl,A,B,s0,s1,l0,l1,theta,t,N);
epMAlt = temp.epM_dephased/(d/(d+1))^M

