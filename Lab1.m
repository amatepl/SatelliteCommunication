clear all; close all; clc;

%Generation of a random bit sequence
n=1e6;
bit_tx = randi(2,1,n)'-1;%stream of bits

%symbol mapping
Nbps = 4; %nb of bits per symbol
modulation = 'qam';
symb_tx= mapping(bit_tx,Nbps,modulation);

%Computation of halfroot cosine filter
Tsym=1e-6;
fsym=1/Tsym; %symbol freq
beta = 0.3; %roll-off factor
M=4; %over-sampling factor
fsamp=M*fsym; %sampling freq
RRCtaps=155; % number of filter taps (odd to have symmetry around the origin)
stepOffset =(1/RRCtaps)*fsamp; % resolution in frequency
highestFreq = stepOffset*(RRCtaps -1)/2; % remove last sample belonging to next period ?
freqGrid = linspace(-highestFreq,highestFreq,RRCtaps);
Hrc=zeros(RRCtaps,1);

for i=1:size(freqGrid,2) 
    f=freqGrid(i);
    if abs(f) <= (1-beta)/(2*Tsym)
        Hrc(i)= Tsym;
    elseif abs(f)<= (1 + beta)/(2*Tsym)
        Hrc(i)=(Tsym/2)*(1+ cos(((pi*Tsym)/beta)*(abs(f)- (1-beta)/(2*Tsym))));
    else 
        Hrc(i)=0;
        
    end
end
Hrrc=sqrt(Hrc); % take the square root of the raised-cosine filter
figure(1)
plot(freqGrid,Hrrc)
xlabel('frequency (Hz)')
ylabel('Amplitude')
title('Frequency response of the RRC filter')

Delta_t = 1/fsamp; %time resolution
t= (-(RRCtaps-1)/2:(RRCtaps-1)/2)*Delta_t;
hrrc=ifft(ifftshift(Hrrc), 'symmetric'); % first need to shift the filter in positive frequencies (Matlab cannot handle negative freq)
hrrc = fftshift(hrrc); % Root-raised cosine filter
hconv=conv(hrrc,hrrc); % convolve the RRC with himself to obtain time domain RC
h2 = hconv/max(hconv); % normalize hrc 
hrrc=hrrc/sqrt(max(hconv)); % scale hrrc accordingly

figure(2)
plot(h2)
xlabel('time (s)')
ylabel('Amplitude')
title('Impulse response of the normalized Raised cosine filter')

figure(3)
plot(t,hrrc)
xlabel('time (s)')
ylabel('Amplitude')
title('Impulse response of the RRC filter')

% TODO : use upsampling and downsampling functions in Matlab
% Convolve the signal with the filter
% Add noise
% convolve at receiver
% remove the first RRCtaps-1 samples and then downsample

% upsampling
upsymb_tx=upsample(symb_tx,M);
% convolution at transmitter with the RRC filter
y=conv(upsymb_tx,hrrc);
% y = y(RRCtaps:end); % discard the first RRCtaps-1 samples
%%% We should not discard taps between the 2 convolutions !!!


%noise
sigE=(trapz(abs(y).^2))*(1/fsamp); %signal energy
Eb=sigE/n; % energy per bit
Eb=Eb/2; % power of a bandpass signal = power of its enveloppe divided by 2
EbN0ratioDB=0:10; %dB
BERmatrix = zeros(length(EbN0ratioDB),1);
BERtool = zeros(length(EbN0ratioDB),1);
for i=1:length(EbN0ratioDB)
    EbN0Ratio = 10^(EbN0ratioDB(i)/10); % convert into non-dB
    N0=Eb/EbN0Ratio;
    NoisePower=2*N0*fsamp;
    noise= sqrt(NoisePower/2)*(randn(1,length(y))+1i*randn(1,length(y)));
    
    % noisy signal
    yn=y+noise';
    
    % convolution with RRC at receiver
    usymb_rx = conv(yn,fliplr(hrrc)); % flip hrrc to obtain matched filter at receiver 
    usymb_rx=usymb_rx(RRCtaps:end-RRCtaps+1); % discard non-causal part (first RRCtaps-1) and additional samples at the end of the signal (last RRCtaps-1)
    symb_rx = downsample(usymb_rx,M); % downsample by M
    bit_rx = demapping(symb_rx,Nbps,modulation); % translate symbols back to bits
    
    
    errVec = bit_tx - bit_rx;
    nbBitsErr = sum(abs(bit_tx-bit_rx))/2; %nb of incorrect bits
    BER= nbBitsErr/n;
    BERmatrix(i)=BER;
    BERtool(i)=berawgn(EbN0Ratio,modulation,2^Nbps);
end

figure(4)
semilogy(EbN0ratioDB,BERmatrix)
title('BER')
xlabel('Eb/N0')
ylabel('BER (log)')
grid on;
hold on;
semilogy(EbN0ratioDB,BERtool)
legend('simulation', 'real curve')

% TODO : do the same for all modulation types

