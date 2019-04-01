function [h] = RRCFilter(beta,Fsymbol,Fsampling, RRCTaps)

% beta - roll-off factor
% Fsymbol - symbol frequency
% Fsampling - sampling frequency
% Nss - Number of symbols
% Nsps - Number of samples per symbol
% RCCTaps - Number of taps
% Delta_t - time resolution (1/Fsampling)
% beta = 0.3;
% Fsymbol = 2e6;
% Fsampling = 4*2e6;
% RRCTaps = 33;

T = 1/Fsymbol;

Delta_t = 1/Fsampling;
t = (-(RRCTaps -1)/2:(RRCTaps -1)/2)*Delta_t;
timesupport = 0.5*(RRCTaps-1)*Delta_t;

stepOffset = (1/RRCTaps)*Fsampling;
highestFreq = stepOffset*(RRCTaps -1)/2;
freqGrid = linspace(-highestFreq,highestFreq,RRCTaps);  % We need an uneven number of samples for matlab


%Fsymbol = 2e6;

% interval_3_1 = freqGrid(abs(freqGrid)> (1+beta)*Fsymbol/2  & freqGrid<0);
% H1 = zeros(size(interval_3_1));
% H = H1;
% interval_2_1 = freqGrid((1-beta)*Fsymbol/2 < abs(freqGrid) & abs(freqGrid)<= (1+beta)*Fsymbol/2 & freqGrid<0);
% H1 = (1 + cos(pi*(abs(interval_2_1)-(1-beta)*Fsymbol/2)/(Fsymbol*beta)))/(Fsymbol*2);
% H = [H,H1];
% interval_1_1 = freqGrid(abs(freqGrid)<= (1-beta)*Fsymbol/2 );
% H1 = ones(size(interval_1_1))*T;
% H = [H,H1];
% %H = H1;
% 
% interval_2_2 = freqGrid((1-beta)*Fsymbol/2 < abs(freqGrid) & abs(freqGrid)<= (1+beta)*Fsymbol/2 & freqGrid>0);
% H1 = (1 + cos(pi*(abs(interval_2_2)-(1-beta)*Fsymbol/2)/(Fsymbol*beta)))/(Fsymbol*2);
% H = [H,H1];
% interval_3_2 = freqGrid(abs(freqGrid)> (1+beta)*Fsymbol/2  & freqGrid>=0);
% H1 = zeros(size(interval_3_2));
% H = [H,H1];
% % Normalization of H
% H = sqrt(H/T);

lcorner = (1 - beta) / (2 * T);
rcorner = (1 + beta) / (2 * T);
H = ones(1, RRCTaps) .* (sqrt(T / 2 * (1 + cos(pi * T / beta * (abs(freqGrid) - lcorner)))));
H(abs(freqGrid) < lcorner) = sqrt(T);
H(abs(freqGrid) > rcorner) = 0;

h = ifft(ifftshift(H),'symmetric');
%h = ifft(fftshift(H));
h = fftshift(h);
h = h/max(sqrt(conv(h,h)));
% size(h);
% figure('Name','h(t)');
% plot(t,h)
% %stem(t,h)
% figure('Name','H(f)');
% plot(freqGrid,H)
