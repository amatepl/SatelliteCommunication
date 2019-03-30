clear all; close all;clc;
Fsymbol = 1e6;
beta = 0.3;
RRCTaps = 33;
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


%SignalEnergy = (trapz(abs(signal).^2))*(1/Fsampling);
%Eb = SignalEnergy/NbEncodedBit;

%--------------------------------
%             |
%             |
%             |
%             V
%-------AWGN----------------------

% %Noise Power
% N0 = Eb/EbN0;
% NoisePower = 2*N0*Fsampling;
% 
% %Noise
% sqrt(NoisePower/2)*(randn(1,length(noise))+ 1i*randn(1,length(noise)));

%signal_noise = noise(signal_tx,Eb_No,Fsampling,Nb);
%--------------------------------
%             |
%             |
%             |
%             V
%---------Root raised cosine filter-----------
%signal_rx = conv(signal_noise,filter);
signal_rx1 = conv(signal_tx,filter);
%size(signal_rx)

%  !!! DISCARD ADDITIONAL SAMPLES AFTER CONVOLUTION !!!
signal_rx = signal_rx1(RRCTaps+M-1:end-RRCTaps+M);
%--------------------------------
%             |
%             |
%             |
%             V
%--------Downsampling------------
down_signal_rx = downsample(signal_rx,M);
%down_signal_rx = signal_rx1;
%size(down_signal_rx)
%--------------------------------
%             |
%             |
%             |
%             V
%---------Decoder----------------
bit_rx = demapping(down_signal_rx,Nbps,'qam');
errors = abs(bit_rx - bit_tx);
sum(errors)
