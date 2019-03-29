
SignalEnergy = (trapz(abs(signal).^2))*(1/Fsampling);
Eb = SignalEnergy/NbEncodedBit;
%-------Encoder------------------
%encoded_signal = mapping(bit_tx,Nbps,modulation);
%--------------------------------
%             |
%             |
%             |
%             V
%-----------Upsampling------------

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


%filter = RRCFilter();
%filtered_signal = conv(upsampled_signal,filter);
%--------------------------------
%             |
%             |
%             |
%             V
%-------AWGN----------------------

%Noise Power
N0 = Eb/EbN0;
NoisePower = 2*N0*Fsampling;

%Noise
sqrt(NoisePower/2)*(randn(1,length(noise))+ 1i*randn(1,length(noise)));
%--------------------------------
%             |
%             |
%             |
%             V
%---------Root raised cosine filter-----------
%filtered_signal = conv(noised_signal,filter);

%  !!! DISCARD ADDITIONAL SAMPLES AFTER CONVOLUTION !!!
%--------------------------------
%             |
%             |
%             |
%             V
%--------Downsampling------------

%--------------------------------
%             |
%             |
%             |
%             V
%---------Decoder----------------

%decoded_signal = demapping(symb_rx,Nbps,modulation);