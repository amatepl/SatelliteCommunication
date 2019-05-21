clear all; close all;clc;
% 1 MHz cutoff fre-quency
% 0.3 roll-off factor  
% Before filtering, the symbol sequence is over-sampled with a factor M >
% 1.
% The sample rate fixes the bandwidth simulated in Matlab, that must be high
% enough to simulate the halfroot Nyquist filtering.
Fsymbol = 2e6;
beta = 0.3;
% !!! WE MUST HAVE ENOUGH TAPS,HOWEVER ERROR OCCURS EVEN WITHOUT NOISE !!!
RRCTaps = 155;
% Upsampling factor 
M = 2;
Fsampling = M*Fsymbol;
Nb = 600000;
% Creation of a random sequence of bits for the moment (row vector)
bit_tx = randi([0 1],Nb, 1);
Nbps = [1 2 4 6];
Eb_No_dB = 1:1:17;
BER = zeros(length(Nbps),length(Eb_No_dB));

% Filter that is used :
filter = RRCFilter(beta,Fsymbol,Fsampling, RRCTaps);
figure();
for j = 1:length(Nbps)
    % TX side :
    [symb_tx,signal_tx] = TX(bit_tx, filter,Nbps(j), M);
    for i = 1:length(Eb_No_dB)
        % Add noise on the signal :
        signal_rx = noise(signal_tx, Eb_No_dB(i),Fsampling,Nb);
        % RX side :
        [symb_rx,bit_rx] = RX(signal_rx, filter,Nbps(j), M, RRCTaps);
        errors = abs(bit_rx - bit_tx);
        BER(j,i) = sum(errors)/Nb;
    end
    semilogy(Eb_No_dB,BER(j,:));
    hold on;
end
title('BER of an ideal channel in satellite communication with AWGN PAM');
legend('BPSK','QPSK','QAM16','QAM64');
xlabel('Eb/N0 (dB)');
ylabel('BER');
grid on;


figure();
berTheo = berawgn(Eb_No_dB,'pam',2^Nbps(1));
semilogy(Eb_No_dB,BER(1,:));
hold on;
semilogy(Eb_No_dB,berTheo);
title('BER of an ideal channel in satellite communication with AWGN PAM');
legend('Simulation','Theorical');
xlabel('Eb/N0 (dB)');
ylabel('BER');
grid on;

figure();
berTheo = berawgn(Eb_No_dB,'qam',2^Nbps(2));
semilogy(Eb_No_dB,BER(2,:));
hold on;
semilogy(Eb_No_dB,berTheo);
title('BER of an ideal channel in satellite communication with AWGN QAM4');
legend('Simulation','Theorical');
xlabel('Eb/N0 (dB)');
ylabel('BER');
grid on;

figure();
berTheo = berawgn(Eb_No_dB,'qam',2^Nbps(3));
semilogy(Eb_No_dB,BER(3,:));
hold on;
semilogy(Eb_No_dB,berTheo);
title('BER of an ideal channel in satellite communication with AWGN QAM16');
legend('Simulation','Theorical');
xlabel('Eb/N0 (dB)');
ylabel('BER');
grid on;

figure();
berTheo = berawgn(Eb_No_dB,'qam',2^Nbps(4));
semilogy(Eb_No_dB,BER(4,:));
hold on;
semilogy(Eb_No_dB,berTheo);
title('BER of an ideal channel in satellite communication with AWGN QAM64');
legend('Simulation','Theorical');
xlabel('Eb/N0 (dB)');
ylabel('BER');
grid on;

