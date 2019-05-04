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
Eb_No_dB = 5;
Nb = Nbps*1e3*blocksize;

% Creation of a random sequence of bits for the moment (row vector)
bit_tx = randi([0 1],1,Nb);

%% LDPC Encoder and generation of the parity check matrix
H = makeLdpc(blocksize,codesize,0,1,3);
% The bits transmitted are made in blocks
bit_tx = reshape(bit_tx,blocksize,Nb/blocksize);
% Generate parity check bits in function of each blocks
[parity_bit,Hnew] = makeParityChk(bit_tx,H,1);
LDPC_bits_tx = vertcat(parity_bit,bit_tx);
LDPC_bits_tx = LDPC_bits_tx(:);

% Filter that is used :
filter = RRCFilter(beta,Fsymbol,Fsampling, RRCTaps);

% TX side :
[symb_tx,signal_tx] = TX(LDPC_bits_tx, filter,Nbps, M);

% Add noise on the signal :
signal_rx = noise(signal_tx, Eb_No_dB,Fsampling,Nb);

% RX side :
[symb_rx,LDPC_bits_rx] = RX(signal_rx, filter,Nbps, M, RRCTaps);
bit_rx = zeros(codesize,Nb/blocksize);
for i = 0:1:Nb/blocksize-1
   bit_rx(:,i+1) = tanner(LDPC_bits_rx(i*codesize+1:(i+1)*codesize), Hnew, 1);
end

bit_rx = bit_rx(blocksize+1:codesize,:);
bit_rx = bit_rx(:);
bit_tx = bit_tx(:);
errors = abs(bit_tx-bit_rx);
BER = sum(errors)/Nb;




