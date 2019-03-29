function [h] = RRCFliter(beta,Fsymbol,Fsampling)

% beta - roll-off factor
% Fsymbol - symbol frequency
% Fsampling - sampling frequency
% Nss - Number of symbols
% Nsps - Number of samples per symbol
% RCCTaps - Number of taps
% Delta_t - time resolution (1/Fsampling)

%Fsymbol = 1e6;
%beta = 0.3;
RRCTaps = 133;
%Fsampling = 8e6;

Delta_t = 1/Fsampling;
t = (-(RRCTaps -1)/2:(RRCTaps -1)/2)*Delta_t;

stepOffset = (1/RRCTaps)*Fsampling;
highestFreq = stepOffset*(RRCTaps -1)/2;
freqGrid = linspace(-highestFreq,highestFreq-1,RRCTaps);  % We need an uneven number of samples for matlab

interval_3_1 = freqGrid(abs(freqGrid)> (1+beta)*Fsymbol/2  & freqGrid<0);
H1 = zeros(size(interval_3_1));
H = H1;
interval_2_1 = freqGrid((1-beta)*Fsymbol/2 < abs(freqGrid) & abs(freqGrid)<= (1+beta)*Fsymbol/2 & freqGrid<0);
H1 = (1 + cos(pi*(abs(interval_2_1)-(1-beta)*Fsymbol/2)/(Fsymbol*beta)))/(Fsymbol*2);
H = [H,H1];
interval_1_1 = freqGrid(abs(freqGrid)<= (1-beta)*Fsymbol/2 );
H1 = ones(size(interval_1_1))/Fsymbol;
H = [H,H1];

interval_2_2 = freqGrid((1-beta)*Fsymbol/2 < abs(freqGrid) & abs(freqGrid)<= (1+beta)*Fsymbol/2 & freqGrid>0);
H1 = (1 + cos(pi*(abs(interval_2_2)-(1-beta)*Fsymbol/2)/(Fsymbol*beta)))/(Fsymbol*2);
H = [H,H1];
interval_3_2 = freqGrid(abs(freqGrid)> (1+beta)*Fsymbol/2  & freqGrid>=0);
H1 = zeros(size(interval_3_2));
H = [H,H1];
H = complex(H);
%H = H + zeros(size(H)).*1i;
H = sqrt(H);

%h = ifft(H,'symmetric');
h = ifft(fftshift(H));
h = ifftshift(h);
%size(h)
figure('Name','real h(t)');
plot(t,h)
%fvtool(h);
h = abs(h)/max(abs(h));
figure('Name','h(t)');
plot(t,h)
figure('Name','H(f)');
plot(freqGrid,H)
grid on;


