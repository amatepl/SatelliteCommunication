function [symb_rx,bit_rx] = RX(signal_rx, filter,Nbps, M, RRCTaps)
% This demodulation uses the demapping function and a downsampling to 
% after a convolution with an ideal filter.
% INPUTS :
% - signal_rx : signal which arrives at the receiver
% - Nbps : number of bits per symbol
% - M : Factor for downsampling
% - filter : row vector in time of the filter which is convolves with
%            the modulated signal. Time domain.
% - RRCTaps : Number of taps
% OUTPUT :
% - symb_rx : complex value of each symbol via the constellation points
% - bit_rx : Interpretation of each symb by encoded bit. 

%% -----------------Filtering-----------------
% Convolution :
up_symb_rx = conv(signal_rx,filter','full');
%up_symb_rx = signal_rx;
% Discard additional samples after convolution
%up_symb_rx = up_symb_rx(RRCTaps + M -1:end-RRCTaps + M);
up_symb_rx = up_symb_rx(RRCTaps:end-RRCTaps+1);
%up_symb_rx = up_symb_rx(:);
%---------------------------------------------
%                      |
%                      |
%                      |
%                      V
%% ---------------Downsampling----------------

symb_rx = downsample(up_symb_rx,M);
%---------------------------------------------
%                      |
%                      |
%                      |
%                      V
%% -----------------Decoder-------------------
if Nbps > 1 
    bit_rx = demapping(symb_rx, Nbps, 'qam');
else 
    bit_rx = demapping(symb_rx, Nbps, 'pam');
end
end
