function signal_rx = noise(signal_tx,Eb_No,Fsampling,Nb)

signal_energy = (trapz(abs(signal_tx).^2))*(1/Fsampling);
Eb = signal_energy/Nb;
Eb = Eb/2;
No = Eb/Eb_No;
noise_power = 2*No*Fsampling;
%noise = ones(length(signal_tx),1);
noise = sqrt(noise_power/2)*(randn(length(signal_tx),1)+1i*randn(length(signal_tx),1));
signal_noise = signal_tx + noise;