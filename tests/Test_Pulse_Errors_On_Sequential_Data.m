%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 29, 2023
%--------------------------------------------------------------------------

load('/Users/evatakou/Documents/MATLAB/GHZ_States_Public/Simulations/GHZ_Data_of_Sequential_FINAL/GHZ4_Sequential.mat')

pulse_err=8/100;

indx=3;

At=OUT.A_Target{indx};
Bt=OUT.B_Target{indx};
t =OUT.Opt_Unit_Times{indx};
N = OUT.Opt_Iters{indx};

%%
[s0,s1,wL]=load_fixed_params;

d=2;
M = length(At)+1;

pref=(d/(d+1))^M;


temp = SuperClass_Sequences(wL,At,Bt,s0,s1,length(At),1,1);
temp = temp.CPMG_error(t,N,pulse_err);
U    = temp.Uevol;


temp=SubClass_Ent_and_Fid;
temp.Uval=U;
temp=temp.MWayEp;



epM_uni = temp.epM_uni/pref

%%
Ze = [s0,0;0,s1];
Z  = [1 0 ; 0 -1];
X  = [0 1 ; 1 0];
Y  = [0 -1i ;1i 0];

IZII=kron(eye(2),kron(Z,eye(4)));
IIZI=kron(eye(4),kron(Z,eye(2)));
IIIZ=kron(eye(8),Z);

ZZII=kron(Ze,kron(Z,eye(4)));
ZIZI=kron(Ze,kron(eye(2),kron(Z,eye(2))));
ZIIZ=kron(Ze,kron(eye(4),Z));

ZXII=kron(Ze,kron(X,eye(4)));
ZIXI=kron(Ze,kron(eye(2),kron(X,eye(2))));
ZIIX=kron(Ze,kron(eye(4),X));

H=wL/2*(IZII+IIZI+IIIZ)+At(1)/2*(ZZII)+At(2)/2*(ZIZI)+At(3)/2*(ZIIZ)+...
                        Bt(1)/2*(ZXII)+Bt(2)/2*(ZIXI)+Bt(3)/2*(ZIIX);

c=2*pi*1e-3;
H=H*c;

Uf=@(t) expm(-1i*H*t);

Rx=@(err) kron(expm(-1i*(pi+err*pi)/2*X),eye(8));

UCPMG=@(t,N,err) (Uf(t/4)*Rx(err)*Uf(t/2)*Rx(err)*Uf(t/4))^N;


UCPMG_tot=1;

for jj=1:length(t)
    
   UCPMG_tot = UCPMG(t(jj),N(jj),pulse_err)*UCPMG_tot ;
    
end


temp2=SubClass_Ent_and_Fid;
temp2.Uval=UCPMG_tot;

temp2=temp2.MWayEp;

ep_M_Alt=temp2.epM_uni;

ep_M_Alt=ep_M_Alt/pref









