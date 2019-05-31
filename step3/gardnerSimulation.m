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
    %signal_rx = noise(signal_tx_M, Eb_No_dB,Fsampling,Nb);
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
    %current = plot(1:step:size(means,2),means(ii,1:step:end));
    current = plot(X,means(:,1:step:end).');
    
    c = get(current, 'color');
%     plot(1:step:size(means,2),means(ii,1:step:end) + stdv(ii,1:step:end),':','Color',c)
%     plot(1:step:size(means,2),means(ii,1:step:end) - stdv(ii,1:step:end),':','Color',c)
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
    
%     
%     symb_rx2 = corrected';
%     bit_rx = zeros(Nb,length(Eb_No_dB));
%     parfor i = 1:length(Eb_No_dB)
%         l_bit_rx = zeros(length(Eb_No_dB),1);
%         if Nbps(j) > 1 
%             l_bit_rx = demapping(symb_rx2(:,i), Nbps(j), 'qam');
%         else 
%             l_bit_rx = demapping(symb_rx2(:,i), 1, 'pam');
%         end 
%         bit_rx(:,i) = l_bit_rx;
%     end
%     
%    
%     
%     
%     
%     
%     errors = abs(bit_rx - kron(bit_tx,ones(1,length(Eb_No_dB))));
%     BER(j,:) = sum(errors)/Nb;
% 
%     semilogy(Eb_No_dB,BER(j,:));
%     hold on;
% %end
% title('Effect of CFO on BER (ISI Only)');
% legend('BPSK','QPSK','QAM16','QAM64');
% xlabel('Eb/N0 (dB)');
% ylabel('BER');
% grid on;
% hold on;
% berTheo1 = berawgn(Eb_No_dB,'qam',2^Nbps(1));
% semilogy(Eb_No_dB,berTheo1,'-o');
% hold on;
% berTheo2 = berawgn(Eb_No_dB,'qam',2^Nbps(2));
% semilogy(Eb_No_dB,berTheo2,'-o');
% hold on;
% berTheo3 = berawgn(Eb_No_dB,'qam',2^Nbps(3));
% semilogy(Eb_No_dB,berTheo3,'-o');
% hold on;
% berTheo4 = berawgn(Eb_No_dB,'qam',2^Nbps(4));
% semilogy(Eb_No_dB,berTheo4,'-o');
% ylim([10e-5 1]);
%legend('BPSK theo','QPSK theo','QAM16 theo','QAM64 theo');

% [symb_tx,signal_tx] = TX(bit_tx, filter,2, M);
% figure;
% plot(symb_tx,'o');
% hold on;
% plot(cfo(signal_rx,df,phi,Fsampling),'x')
% legend('Emitted symbols', 'phase offset');
% title('Effect of phase offset pi/16')
% grid on;


    