clear all; close all; clc;
Fsymbol = 1e6;
beta = 0.3;
RRCTaps = 133;
Fsampling = 8e6;
Nb = 10000;

% Creation of a random sequence of bits for the moment (row vector)
bit_tx = randi([0 1],Nb, 1);
Nbps = 2; % Bit per symbol 
% Upsampling factor 
M = 4;

Eb_No = 2;
filter = RRCFilter(beta,Fsymbol,Fsampling, RRCTaps);
signal_tx = TX(bit_tx,filter,Nbps,M);
signal_noise = noise(signal_tx,Eb_No,Fsampling,Nb);
%size(signal_noise);
