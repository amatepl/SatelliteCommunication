function [symb_tx,signal_tx_kron,signal_tx_up] = TX(bit_tx, filter,Nbps, M)
% This modulation uses the mapping function and an oversampling to 
% make the convolution with an ideal filter.
% INPUTS :
% - bit_tx : vector of input bits 
% - Nbps : number of bits per symbol
% - M : Factor for oversampling
% - filter : row vector in time of the filter which is convolves with
%            the modulated signal. Time domain.
% OUTPUT :
% - signal_TX : Signal that is sended by the transmitter

%% ------------------Encoder------------------
% Digital modulation (use qam or pam => Depend of the Nbps)
if Nbps > 1 
    symb_tx = mapping(bit_tx, Nbps, 'qam');
else 
    symb_tx = mapping(bit_tx, Nbps, 'pam');
end
%---------------------------------------------
%                       |
%                       |
%                       |
%                       V
%% -----------------Upsampling----------------
% Upsampling vector :
% M times the same sample of symb_tx
up_symb_tx_kron = kron(symb_tx,ones([M 1]));
up_symb_tx_up = upsample(symb_tx,M);
%---------------------------------------------
%                       |
%                       |
%                       |
%                       V
%% -----------------Filtering-----------------
% Convolution :
signal_tx_kron = conv(up_symb_tx_kron,filter);
signal_tx_up = conv(up_symb_tx_up,filter);
end


