function [signal_rx,No] = noise(signal_tx,Eb_No_dB,Fsampling,Nb)
% Add an additive Gaussian White noise (AWGN) on the signal inside the
% channel.
% INPUTS :
% - signal_tx : signal which is transmitted
% - Eb_No : Energy ratio between density energy of one bit and the density
%           energy of the noise.
% - Fsampling : sampling frequency
% - Nb : Number of bits
% OUTPUTS :
% - signal_rx : signal which arrives at the receiver


Eb_No = 10.^(Eb_No_dB/10);
signal_energy = (trapz(abs(signal_tx).^2))*(1/Fsampling);
Eb = signal_energy/Nb;
Eb = Eb/2;
No = Eb/Eb_No;
noise_power = 2*No*Fsampling;
%noise_power = kron(noise_power,ones(length(signal_tx),1));

noise = randn(length(signal_tx),1)+1i*randn(length(signal_tx),1);
%noise = kron(noise,ones([1,length(No)]));

noise = sqrt(noise_power/2)*noise;
signal_rx = signal_tx + noise;