
function signal_tx = TX(bit_tx, filter,Nbps, M)
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

% Digital modulation (use qam or pam => Depend of the Nbps)
if Nbps > 1 
    symb_tx = mapping(bit_tx, Nbps, 'qam');
else 
    symb_tx = mapping(bit_tx, Nbps, 'pam');
end

% Upsampling vector :
up_symb_tx = upsample(symb_tx,M);

% Convolution :
signal_tx = conv(up_symb_tx,filter);
end

