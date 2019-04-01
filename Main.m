clear all; close all;clc;
Fsymbol = 2e6;
beta = 0.3;
% !!! WE MUST HAVE ENOUGH TAPS,HOWEVER ERROR OCCURS EVEN WITHOUT NOISE !!!
RRCTaps = 131;

% Upsampling factor 
M = 4;

Fsampling = M*2e6;

Nb = 120000;

% Creation of a random sequence of bits for the moment (row vector)
bit_tx = randi([0 1],Nb, 1);
Eb_No_dB = 1:1:50;

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
%filter = RRCFilter(beta,Fsymbol,Fsampling, RRCTaps).';
filter = RRCFilter(beta,Fsymbol,Fsampling, RRCTaps);



for Nbps = [1 2 4 6]
    % TX side :
    [symb_tx,signal_tx_up] = TX(bit_tx, filter,Nbps, M);
    BER_up = zeros(1,length(Eb_No_dB));
    for i = 1:length(Eb_No_dB)
        % Add noise on the signal :
        signal_rx_up = noise(signal_tx_up, Eb_No_dB(i),Fsampling,Nb);
%                     |
%                     |
%                     V
        % RX side :
        [symb_rx_up,bit_rx_up] = RX(signal_rx_up, filter,Nbps, M, RRCTaps);
        errors_up = abs(bit_rx_up - bit_tx);
        BER_up(i) = sum(errors_up)/2/Nb;
    end
    figure(3);
    semilogy(Eb_No_dB, BER_up);
    hold on;
    title('BER of an ideal channel in satellite communication with AWGN');
    legend('BPSK','QPSK','16QUAM','64QUAM');
    xlabel('Eb/N0 (dB)');
    ylabel('BER');
end
