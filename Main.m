clear all; close all;clc;
Fsymbol = 1e6;
beta = 0.3;
% !!! WE MUST HAVE ENOUGH TAPS,HOWEVER ERROR OCCURS EVEN WITHOUT NOISE !!!
RRCTaps = 133;
Fsampling = 4e6;
Nb = 200000;

% Creation of a random sequence of bits for the moment (row vector)
bit_tx = randi([0 1],Nb, 1);
Nbps = 2; % Bit per symbol 
% Upsampling factor 
M = Fsampling/Fsymbol;

Eb_No_dB = [1:1:50];
BER_kron = zeros(1,length(Eb_No_dB));
BER_up = zeros(1,length(Eb_No_dB));

% 1 MHz cutoff fre-quency
% 0.3 roll-off factor  
% Before filtering, the symbol sequence is over-sampled with a factor M >
% 1.
% The sample rate fixes the bandwidth simulated in Matlab, that must be high
% enough to simulate the halfroot Nyquist filtering.

% Filter that is used :
filter = RRCFilter(beta,Fsymbol,Fsampling, RRCTaps).';
% TX side :
[symb_tx,signal_tx_kron,signal_tx_up] = TX(bit_tx, filter,Nbps, M);
for i = 1:length(Eb_No_dB)
    % Add noise on the signal :
    signal_rx_kron = noise(signal_tx_kron, Eb_No_dB(i),Fsampling,Nb);
    signal_rx_up = noise(signal_tx_up, Eb_No_dB(i),Fsampling,Nb);
    % signal_rx_kron = signal_tx;
    % signal_rx_up = signal_tx;
    % RX side :
    [symb_rx_kron,bit_rx_kron] = RX(signal_rx_kron, filter,Nbps, M, RRCTaps);
    [symb_rx_up,bit_rx_up] = RX(signal_rx_up, filter,Nbps, M, RRCTaps);
    errors_kron = abs(bit_rx_kron - bit_tx);
    errors_up   = abs(bit_rx_up - bit_tx);
    BER_kron(i) = sum(errors_kron)/2/Nb;
    BER_up(i) = sum(errors_up)/2/Nb;
end

figure(3);
semilogy(Eb_No_dB, BER_kron);
hold on;
semilogy(Eb_No_dB, BER_up);
title('BER of an ideal channel in satellite communication with AWGN');
legend('kronecker','upsample');
xlabel('Eb/N0 (dB)');
ylabel('BER');

figure(4);
plot(symb_rx_up, 'o');
title('symbol at receiver with upsample');

figure(5);
plot(symb_rx_kron, 'o');
title('symbol at receiver with kronecker');

