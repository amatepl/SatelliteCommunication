clear all; close all; clc;
Fsymbol = 1e6;
beta = 0.3;
RRCTaps = 133;
Fsampling = 8e6;

% Creation of a random sequence of bits for the moment (row vector)
bit_tx = randi([0 1],10000, 1);
Nbps = 2; % Bit per symbol 
% Upsampling factor 
M = 4;
filter = RRCFilter(beta,Fsymbol,Fsampling, RRCTaps);
signal_tx = TX(bit_tx,filter,Nbps,M);