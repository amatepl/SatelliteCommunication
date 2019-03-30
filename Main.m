clear all; close all;clc;
Fsymbol = 1e6;
beta = 0.3;
% !!! WE MUST HAVE ENOUGH TAPS,HOWEVER ERROR OCCURS EVEN WITHOUT NOISE !!!
RRCTaps = 133;
Fsampling = 4e6;
Nb = 10000;

% Creation of a random sequence of bits for the moment (row vector)
bit_tx = randi([0 1],Nb, 1);
Nbps = 2; % Bit per symbol 
% Upsampling factor 
M = Fsampling/Fsymbol;

Eb_No = 2;

% 1 MHz cutoff fre-quency
% 0.3 roll-off factor  
% Before filtering, the symbol sequence is over-sampled with a factor M >
% 1.
% The sample rate fixes the bandwidth simulated in Matlab, that must be high
% enough to simulate the halfroot Nyquist filtering.

% Filter that is used :
filter = RRCFilter(beta,Fsymbol,Fsampling, RRCTaps).';
% TX side :
[symb_tx,signal_tx] = TX(bit_tx, filter,Nbps, M);
% Add noise on the signal :
% signal_rx = noise(signal_tx, Eb_No,Fsampling,Nb);
signal_rx = signal_tx;
% RX side :
[symb_rx,bit_rx] = RX(signal_rx, filter,Nbps, M, RRCTaps);

errors = abs(bit_rx - bit_tx);
sum(errors)
