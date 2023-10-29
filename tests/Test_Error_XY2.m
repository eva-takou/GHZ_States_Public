%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 28, 2023
%%
clc;
clear;
close all

%% Define parameters of a single spin and get the evolution operator based on the class
A=40;
B=60;
wL=314;
s0=0;
s1=-1;



t=9.8;
N=6;

pulse_error=0.08;



spin=SuperClass_Sequences(wL,A,B,s0,s1,1,1,1);
spin=spin.XY2_error(t,N,pulse_error);


U=spin.Uevol;

%% Get the evolution operator in the following way:

c=2*pi*1e-3;


X=[0 1 ; 1 0];
Y=[0 -1i ; 1i 0];
Ze=[s0 0 ; 0 s1];
Z=[1 0 ; 0 -1];

IZ=kron(eye(2),Z);
ZZ=kron(Ze,Z);
ZX=kron(Ze,X);


Rx=@(phi) kron(expm(-1i*phi/2*X),eye(2));
Ry=@(phi) kron(expm(-1i*phi/2*Y),eye(2));

H = (wL/2*IZ+A/2*ZZ+B/2*ZX)*c;

Ufree=@(t) expm(-1i*t*H);

UXY2 =@(err) Ufree(t/4)*Ry(pi+err*pi)*Ufree(t/2)*Rx(pi+err*pi)*Ufree(t/4);

Ualt = UXY2(pulse_error)^N;

Ualt-U










