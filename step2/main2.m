addpath('../shared');
addpath('../step1');
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
Nbps = 1;
Eb_No_dB = 1:1:7;
BER = zeros(length(Nbps),length(Eb_No_dB));

% LDPC parameters :
blocksize = 126;
coderate = 1/2;
codesize = blocksize/coderate;
Nb = Nbps*1e2*blocksize;
maxit = [1,5,10];

for j = 1:length(maxit)
    % Parity check matrix
    H = makeLdpc(blocksize,codesize,0,1,4);
    block_bit = randi([0 1],Nb,1);
    save = block_bit;
    block_bit = reshape(block_bit,blocksize,Nb/blocksize);
    [parity_bit,Hnew] = makeParityChk(block_bit,H,1);
    bit_tx = [parity_bit;block_bit];
    bit_tx = bit_tx(:);
    % Filter that is used :
    filter = RRCFilter(beta,Fsymbol,Fsampling, RRCTaps);
    % TX side :
    [symb_tx,signal_tx] = TX(bit_tx, filter,Nbps, M);
    tic;
    for i = 1:length(Eb_No_dB)
        % Add noise on the signal :
        [signal_rx,No] = noise(signal_tx, Eb_No_dB(i),Fsampling,length(bit_tx));
        % RX side :
        [symb_rx,bit_rx] = RX(signal_rx, filter,Nbps, M, RRCTaps);
%         for k = 0:codesize:length(bit_rx)-codesize
%             bit_rx(k+1:k+codesize) = LDPC(Hnew,bit_rx(k+1:k+codesize).',maxit(j));
%         end
        errors = abs(bit_rx - bit_tx);    
        BER(j,i) = sum(errors)/length(errors);
    end
    toc;
end

figure();
berTheo = berawgn(Eb_No_dB,'pam',2^Nbps);
for j = 1:length(maxit)
    semilogy(Eb_No_dB,BER(j,:),'-o');
    hold on;
end
semilogy(Eb_No_dB,berTheo,'-o');
title('BER for hard decoding of PAM modulation');
legend('Simulation 1 it','Simulation 5 it','Simulation 10 it','Theorical');
xlabel('Eb/N0 (dB)');
ylabel('BER');
grid on;

