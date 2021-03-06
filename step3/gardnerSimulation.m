  
clear all; clc; close all; clear global;
addpath('../step1');
addpath('../shared');

% Channel parameters
Fsymbol = 1e6;                          % Symbol frequency
M = 30;                                 % Uppsampling factor
Fsampling = M*Fsymbol;                  % Sampling frequency
RRCTaps = 151;                          % RCC filter taps
beta = 0.3;                             % Roll-off factor
Fgardner = Fsampling/2;                 % Gardner frequency
fc=2e9;                                 % Carrier frequency

% Simualtion parameters----------------------------------------------------

%   Dose are the parameters to play with.
k=[0.05 0.1 0.4];                                % kappa - correction weight factor
%k=[0.1];
shift = round(0.45*Fsampling/Fsymbol);  % time shift
Nbps = [6];                             % Number of bits per symbol
Nb = 2000*Nbps;                         % Number of bits
Eb_No_dB = 0:2:16;                      % Bit energy over noise energy in dB
%df= (1e-6*fc)*[0 10 50];                     % Carrier frequency offset
df = (1e-6*fc)*[10];
phi=0;                                  % Phase offset

%--------------------------------------------------------------------------

%copies = ones(1,length(Eb_No_dB));
%copies = ones(1,length(df));
copies = ones(1,length(k));
%BER = zeros(length(Nbps),length(Eb_No_dB));

% Nyquist
filter = RRCFilter(beta,Fsymbol,Fsampling, RRCTaps);
%figure;
j=1;
%for j = 1:length(Nbps)
Nb_symb = Nb/Nbps;
means = zeros(length(k), Nb_symb);
stdv = zeros(length(k), Nb_symb);
step = 20;
figure('Name','Gardner')
hold on
ii =1;
Nbexp = 50;
for jj = 1:Nbexp
    
    % Creation of a random sequence of bits for the moment (row vector)
    bit_tx = randi([0 1],Nb, 1);
    
    % TX side :
    [symb_tx,signal_tx] = TX(bit_tx, filter,Nbps(j), M);
    signal_tx_M =  kron(signal_tx,copies);

    % Add noise on the signal :
    %signal_rx = noise2(signal_tx_M, Eb_No_dB,Fsampling,Nb);
    signal_rx = signal_tx_M;
    %num = size(signal_rx,1);
    % Add CFO
    t = (0:length(signal_rx)-1)./Fsampling;
    
    %Use this when plotting for Eb/No
    %signal_cfo =signal_rx .* kron(exp(1j.* (2*pi*df .* t' + phi*ones(length(signal_rx),1))),copies);
    
    %Use this when plotting for CFO
    signal_cfo =signal_rx .* exp(1j.* (2*pi*df .* t' + phi*ones(length(signal_rx),1)));
    %clear signal_rx;

    % RX side :
    symb = zeros(Nb*M/Nbps,size(signal_cfo,2));
   
    parfor n = 1:size(signal_cfo,2)
        up_symb_rx = conv(signal_cfo(:,n),filter','full');
        symb(:,n) = up_symb_rx(RRCTaps:end-RRCTaps+1);
    end
    symb = symb(1+shift:end,:);

    size_symb = size(symb);
    t = (0:size_symb(1)-1)./Fsampling;
    
    %clear signal_cfo;
    symb_rx = symb;
    
    [~, epsilon] = gardner(symb_rx,k,M);
    samplingerror = (shift - epsilon*Fsampling/Fsymbol)/Fsampling;
    means = means + samplingerror;
    stdv = stdv + samplingerror.^2;
end
    
    means = means/Nbexp;
    stdv = sqrt(stdv/Nbexp - means.^2);
    X = kron([1:step:size(means,2)].',copies);
    current = plot(X,means(:,1:step:end).');
    
    c = get(current, 'color');
    plot_stdv_up = plot(X,means(:,1:step:end).' + stdv(:,1:step:end).',':','LineWidth',2);
    set(plot_stdv_up,{'color'},c);
    plot_stdv_down = plot(X,means(:,1:step:end).' - stdv(:,1:step:end).',':','LineWidth',2);    
    set(plot_stdv_down,{'color'},c);
%     stdv =  sum(stdv,2)/1000;
%     plot_noise = plot(Eb_No_dB,stdv.');
    legend(current,{'K = 0.01','K = 0.1','K = 0.4'});
    xlabel('Symbols');
    ylabel('Time error (mean +- stdv)');
    grid on;
    