function [h] = RRCFilter(beta,Fsymbol,Fsampling, RRCTaps)
% Root raised cosine filter :
% INPUTS :
% - beta : roll-off factor
% - Fsymbol : symbol frequency
% - Fsampling : sampling frequency
% - RCCTaps : Number of taps
% OUTPUTS :
% - h : the normalized root raised cosine filter in the time domain
T = 1/Fsymbol;

Delta_t = 1/Fsampling;
t = (-(RRCTaps -1)/2:(RRCTaps -1)/2)*Delta_t;
timesupport = 0.5*(RRCTaps-1)*Delta_t;

stepOffset = (1/RRCTaps)*Fsampling;
highestFreq = stepOffset*(RRCTaps -1)/2;
freqGrid = linspace(-highestFreq,highestFreq,RRCTaps);  % We need an uneven number of samples for matlab

interval_3_1 = freqGrid(abs(freqGrid)> (1+beta)*Fsymbol/2  & freqGrid<0);
H1 = zeros(size(interval_3_1));
H = H1;
interval_2_1 = freqGrid((1-beta)*Fsymbol/2 < abs(freqGrid) & abs(freqGrid)<= (1+beta)*Fsymbol/2 & freqGrid<0);
H1 = (1 + cos(pi*(abs(interval_2_1)-(1-beta)*Fsymbol/2)/(Fsymbol*beta)))/(Fsymbol*2);
H = [H,H1];
interval_1_1 = freqGrid(abs(freqGrid)<= (1-beta)*Fsymbol/2 );
H1 = ones(size(interval_1_1))*T;
H = [H,H1];

interval_2_2 = freqGrid((1-beta)*Fsymbol/2 < abs(freqGrid) & abs(freqGrid)<= (1+beta)*Fsymbol/2 & freqGrid>0);
H1 = (1 + cos(pi*(abs(interval_2_2)-(1-beta)*Fsymbol/2)/(Fsymbol*beta)))/(Fsymbol*2);
H = [H,H1];
interval_3_2 = freqGrid(abs(freqGrid)> (1+beta)*Fsymbol/2  & freqGrid>=0);
H1 = zeros(size(interval_3_2));
H = [H,H1];

H = sqrt(H);

h = ifft(ifftshift(H),'symmetric');
h = fftshift(h);
h = h/max(sqrt(conv(h,h)));

figure();
plot(freqGrid,H)
xlabel('frequency (Hz)')
ylabel('Amplitude')
title('Frequency response of the RRC filter')

figure();
plot(t,h)
xlabel('time (s)')
ylabel('Amplitude')
title('Time response of the RRC filter with normalization')