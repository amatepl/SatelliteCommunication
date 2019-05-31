function signal_rx = noise2(signal_tx,Eb_No_dB,Fsampling,Nb)
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
signal_energy = (trapz(abs(signal_tx(:,1)).^2))*(1/Fsampling);
Eb = signal_energy/Nb;
Eb = Eb/2;
No = Eb./Eb_No;
size_signal = size(signal_tx);

% Matrix with all the columns = signal_tx and number of columns = length(Eb_No_dB)
clone = ones(size_signal(1),1);
noise_power = kron(2*No*Fsampling,clone);       % lines(signla)xlength(Eb_No_dB)

%noise_power = kron(noise_power,ones(length(signal_tx),1));

%noise = randn(length(signal_tx),1)+1i*randn(length(signal_tx),1);
noise = randn(size(noise_power))+1i*randn(size(noise_power));
%noise = kron(noise,ones([1,length(No)]));

noise = sqrt(noise_power./2).*noise;

signal_rx = signal_tx + noise;