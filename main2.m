clear all; close all;clc;
%% Definition of the variables
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

% we want divide our bit frame into equivalent sub blocks
blocksize = 128; % = number of rows in the parity check matrix
% We want a code rate smaller than one 
coderate = 1/2;  % ex : 1/2 => 1 check bit for 1 bit
                 %      3/4 => 1 check bit for 3 bits
                 % if code rate increase => protection against errors
                 % decreases
                 
% bits transmitted per block = number of columns in the parity check matrix
codesize = blocksize/coderate;

Nbps = 1;
Eb_No_dB = 2:1:5;
Nb = Nbps*1e3*blocksize;

%% LDPC Encoder and generation of the parity check matrix
H = makeLdpc(blocksize,codesize,0,1,3);
% The bits transmitted are made in blocks
% Creation of a random sequence of bits for the moment (row vector)
bit_tx = randi([0 1],blocksize,Nb/blocksize);
% Generate parity check bits in function of each blocks
[parity_bit,Hnew] = makeParityChk(bit_tx,H,1);
LDPC_bits_tx = [parity_bit;bit_tx];

%% Filter + noise :
% Filter that is used :
filter = RRCFilter(beta,Fsymbol,Fsampling, RRCTaps);

BER = zeros(1,length(Eb_No_dB));
bits_tx = LDPC_bits_tx(:);

for j = 1:length(Eb_No_dB)
    LDPC_bits_rx = zeros(codesize,Nb/blocksize);
    for k = 1:length(LDPC_bits_rx)
        % TX side :
        [symb_tx,signal_tx] = TX(LDPC_bits_tx(:,k), filter,Nbps, M);

        % Add noise on the signal :
        signal_rx = noise(signal_tx, Eb_No_dB(j),Fsampling,codesize);

        % RX side :
        [symb_rx,bits_rx] = RX(signal_rx, filter,Nbps, M, RRCTaps);
        LDPC_bits_rx(:,k) = tanner(bits_rx,Hnew,2);
    end
    bits_rx = LDPC_bits_rx(:);
    error= abs(bits_tx-bits_rx);
    BER(j) = 1e6*(sum(error)*coderate)*1/Nb*1e-6;
end
figure();
berTheo = berawgn(Eb_No_dB,'pam',2);
semilogy(Eb_No_dB,BER_coded,'-o');
hold on;
% semilogy(Eb_No_dB,BER_uncoded,'-o');
% hold on;
semilogy(Eb_No_dB,berTheo,'-o');
title('BER of an ideal channel in satellite communication with AWGN PAM');
legend('LDPC hard decoding coded 20 it','Theorical','Location', 'Best');
xlabel('Eb/N0 (dB)');
ylabel('BER');
grid on;


