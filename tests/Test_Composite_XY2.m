%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 27, 2023
%% Evaluate a composite CPMG evolution operator for 1 nuclear spin
clc;
clear;
close all;

t1=5;
t2=12;

t=[t1,t2];
N1=100;
N2=30;
N=[N1,N2];
Nnuc=1;
A=30;
B=60;

k=1;
wL=314;
s0=0;
s1=-1;



spin=SubClass_U4Operations(wL,A,B,s0,s1,Nnuc,k,1);
spin=spin.XY2(t,N);
U=spin.Uevol;

%% Check that the evolution agrees with 2 compositions

spin=SubClass_U4Operations(wL,A,B,s0,s1,Nnuc,k,N1);
spin=spin.XY2(t1);
U1=spin.Uevol;

spin=SubClass_U4Operations(wL,A,B,s0,s1,Nnuc,k,N2);
spin=spin.XY2(t2);
U2=spin.Uevol;

all(all(U2*U1==U))




