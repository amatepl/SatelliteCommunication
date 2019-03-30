clear all; close all;clc;
Fsymbol = 2e6;
beta = 0.3;
% !!! WE MUST HAVE ENOUGH TAPS,HOWEVER ERROR OCCURS EVEN WITHOUT NOISE !!!
RRCTaps = 133;

% Upsampling factor 
M = 4;

Fsampling = M*2e6;

Nb = 240000;

% Creation of a random sequence of bits for the moment (row vector)
bit_tx = randi([0 1],Nb, 1);
Nbps = 2; % Bit per symbol 

Eb_No_dB = 1:1:50;
BER_up = zeros(1,length(Eb_No_dB));

%SignalEnergy = (trapz(abs(signal).^2))*(1/Fsampling);
%Eb = SignalEnergy/NbEncodedBit;

%-------Encoder------------------
%symb_tx = mapping(bit_tx, Nbps, 'qam');
%size(symb_tx)
%--------------------------------
%             |
%             |
%             |
%             V
%-----------Upsampling------------


%up_symb_tx = upsample(symb_tx,M);

%up_symb_tx = kron(symb_tx,ones([M 1]));
%size(up_symb_tx);
%--------------------------------
%             |
%             |
%             |
%             V

%-------Nyquist Filter------------

% 1 MHz cutoff fre-quency
% 0.3 roll-off factor  
% Before filtering, the symbol sequence is over-sampled with a factor M >
% 1.
% The sample rate fixes the bandwidth simulated in Matlab, that must be high
% enough to simulate the halfroot Nyquist filtering.

% Filter that is used :
filter = RRCFilter(beta,Fsymbol,Fsampling, RRCTaps).';


% TX side :



%Eb_No_dB = 1:50/length(signal_tx_up):50;

for Nbps = [1 2 4 6]
    [symb_tx,signal_tx_up] = TX(bit_tx, filter,Nbps, M);
    BER_up = zeros(1,length(Eb_No_dB));
    for i = 1:length(Eb_No_dB)
        % Add noise on the signal :
        signal_rx_up = noise(signal_tx_up, Eb_No_dB(i),Fsampling,Nb);

        % signal_rx = signal_tx;
        % RX side :

        [symb_rx_up,bit_rx_up] = RX(signal_rx_up, filter,Nbps, M, RRCTaps);
        errors_up   = abs(bit_rx_up - bit_tx);
        BER_up(i) = sum(errors_up)/2/Nb;
    end
    figure(3);
    % semilogy(Eb_No_dB, BER_kron);
    % hold on;
    semilogy(Eb_No_dB, BER_up);
    hold on;
    title('BER of an ideal channel in satellite communication with AWGN');
    legend('BPSK','QPSK','16QUAM','64QUAM');
    xlabel('Eb/N0 (dB)');
    ylabel('BER');
end
% figure(3);
% % semilogy(Eb_No_dB, BER_kron);
% % hold on;
% semilogy(Eb_No_dB, BER_up);
% title('BER of an ideal channel in satellite communication with AWGN');
% legend('upsample');
% xlabel('Eb/N0 (dB)');
% ylabel('BER');
