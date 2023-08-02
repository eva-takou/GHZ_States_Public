function [s0,s1,wL]=load_fixed_params
%Script to load the fixed parameters for the simulations in the GHZ paper.

%=== Fixed parameters ===================================================
s0 = 0;  s1 = -1;                   % Spin projections of NV

B0        = 403;                    % External B-field in Gauss
B0        = B0*1e-4;                % B-field in T

gamma_C13 = 6.728284  * 1e7;        % Gyromagnetic ratio rad/T*s
wL        = gamma_C13 * B0 ;        % Larmor freq of nuclei
wL        = wL*1e-6;                % ~2.7 MHz
wL        = wL*1e3;                 % in KHz
wL        = wL/(2*pi);              % 431.5484 kHz


end