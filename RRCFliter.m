function [filtered] = RRCFliter(beta,Ns,Nsps)

% beta - roll-off factor
% Fsymbol - symbol frequency
% Fsampling - sampling frequency
% Nss - Number of symbols
% Nsps - Number of samples per symbol
% RCCTaps - Number of taps
% Delta_t - time resolution (1/Fsampling)


stepOffset = (1/RRCTaps)*fs;
highestFreq = stepOffset*(RCCTaps -1)/2;
freqGrid = linspace(-highestFreq,highestFreq,RCCTaps);

H = 1/(Fsymbol*2)*(1 + cos(pi*(abs(f)-(1-beta)*Fsymbol/2)/(Fsymbol*beta)));

