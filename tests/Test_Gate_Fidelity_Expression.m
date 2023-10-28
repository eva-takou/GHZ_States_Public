clear;
clc;
close all;

Nnuc=5;

A=10+(200-10)*rand(1,Nnuc);
B=10+(200-10)*rand(1,Nnuc);

Target   = [1,2];

Ntarget  = length(Target);
Unwanted = setxor(1:Nnuc,Target);

At=A(Target);
Bt=B(Target);

Au=A(Unwanted);
Bu=B(Unwanted);

wL=314;
s0=0;
s1=-1;

k=1;
N=10;
t=5.6;

for jj=1:length(Au)
    
    spin=SubClass_U4Operations(wL,Au(jj),Bu(jj),s0,s1,1,1,N);
    spin=spin.CPMG(t);
    spin=spin.Rot_Angles_And_Axes;
    n0u(jj,:)=spin.axes{1};
    n1u(jj,:)=spin.axes{2};
    phi0u(jj)=spin.angles{1};
    phi1u(jj)=spin.angles{2};
    
end

temp=SubClass_Ent_and_Fid;
Infid1=temp.Gate_Infid(Nnuc,Ntarget,phi0u,phi1u,n0u,n1u).Infid;

%% Get the target evolution

temp=SuperClass_Sequences(wL,At,Bt,s0,s1,length(At),1,N);
temp=temp.CPMG(t);
U0=temp.Uevol;



%% Now calculate full evolution

temp=SuperClass_Sequences(wL,A,B,s0,s1,Nnuc,1,N);
temp=temp.CPMG(t);
U=temp.Uevol;

temp=SubClass_Ent_and_Fid;

Ek=temp.Get_Kraus(U,Nnuc,Target);


Infid_Alt=temp.Get_Infid_From_Kraus(U0,Ek,Target)


abs(Infid_Alt-Infid1)







 