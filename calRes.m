%%
clear;
close all;

% Inductor parameters
indMeta = zeros(7,1);
indMeta(1) = 6e-3; % OD
indMeta(2) = 0.4e-3; % LW
indMeta(3) = 0.4e-3; % LS
indMeta(4) = 3; % N
indMeta(5) = 0.6e-3; % SEP
indMeta(6) = 0.05e-3; %dl
indMeta(7) = 3.9; % epsR
% permittivity of surrounding media is set in calC_mat.m

L = calL(indMeta);
C = calC_mat(indMeta);

omega0 = 1/sqrt(L*C);
f0 = omega0/(2*pi)