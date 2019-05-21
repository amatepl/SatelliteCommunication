clear; clc; close all;
addpath('../step1');
addpath('../');
% Simulation parameters
Fsymbol = 2e6;
M = 2;
Fsampling = M*Fsymbol;
RRCTaps = 155;
beta = 0.3;
Nbps = [1 2 4 6];

fc=2e9;

Nb = 600000;
% Creation of a random sequence of bits for the moment (row vector)
bit_tx = randi([0 1],Nb, 1);

Eb_No_dB = 1:1:16;
BER = zeros(length(Nbps),length(Eb_No_dB));

df = (1e-6*fc)*25; % ppm
phi=pi/16;
% Nyquist
filter = RRCFilter(beta,Fsymbol,Fsampling, RRCTaps);
figure;
for j = 1:length(Nbps)
    % TX side :
    [~,signal_tx] = TX(bit_tx, filter,Nbps(j), M);
    for i = 1:length(Eb_No_dB)
        % Add noise on the signal :
        signal_rx = noise(signal_tx, Eb_No_dB(i),Fsampling,Nb);

        % Add CFO
        t = (0:length(signal_rx)-1)./Fsampling;
        signal_cfo =signal_rx .* exp(1j.* (2*pi*df .* t'));
        
        % RX side :
        up_symb_rx = conv(signal_cfo,filter','full');
        symb = up_symb_rx(RRCTaps:end-RRCTaps+1);
        t = (0:length(symb)-1)./Fsampling;
        symb = (symb .* exp(1j.* (-2*pi*df .* t'))).*exp(-1j*(RRCTaps-1)/2/Fsampling*2*pi*df);
        
        symb_rx = downsample(symb,M);

        if Nbps(j) > 1 
            bit_rx = demapping(symb_rx, Nbps(j), 'qam');
        else 
            bit_rx = demapping(symb_rx, 1, 'pam');
        end  

        errors = abs(bit_rx - bit_tx);
        BER(j,i) = sum(errors)/Nb;
    end
    semilogy(Eb_No_dB,BER(j,:));
    hold on;
end
title('Effect of CFO on BER (ISI Only)');
legend('BPSK','QPSK','QAM16','QAM64');
xlabel('Eb/N0 (dB)');
ylabel('BER');
grid on;
hold on;
berTheo1 = berawgn(Eb_No_dB,'pam',2^Nbps(1));
semilogy(Eb_No_dB,berTheo1,'-o');
hold on;
berTheo2 = berawgn(Eb_No_dB,'qam',2^Nbps(2));
semilogy(Eb_No_dB,berTheo2,'-o');
hold on;
berTheo3 = berawgn(Eb_No_dB,'qam',2^Nbps(3));
semilogy(Eb_No_dB,berTheo3,'-o');
hold on;
berTheo4 = berawgn(Eb_No_dB,'qam',2^Nbps(4));
semilogy(Eb_No_dB,berTheo4,'-o');
ylim([10e-5 1]);

% [symb_tx,signal_tx] = TX(bit_tx, filter,2, M);
% figure;
% plot(symb_tx,'o');
% hold on;
% plot(cfo(signal_rx,df,phi,Fsampling),'x')
% legend('Emitted symbols', 'phase offset');
% title('Effect of phase offset pi/16')
% grid on;


    