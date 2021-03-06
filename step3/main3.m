clear all; clc; close all; clear global;
addpath('../step1');
addpath('../shared');

% Simulation parameters
Fsymbol = 1e6;
M = 16;
Fsampling = M*Fsymbol;
RRCTaps = 151;
beta = 0.3;
Nbps = [6];
Fgardner = Fsampling/2;
k=0.1;
shift = round(0.45*Fsampling/Fsymbol);
fc=2e9;

Nb = 1000*Nbps;


Eb_No_dB = 4:4:16;
copies = ones(1,length(Eb_No_dB));
BER = zeros(length(Nbps),length(Eb_No_dB));

%df = (1e-6*fc)*25; % ppm
df= (1e-6*fc)*10; % ppm
%df = 0; % ppm
phi=0;
% Nyquist
filter = RRCFilter(beta,Fsymbol,Fsampling, RRCTaps);
%figure;
j=1;
%for j = 1:length(Nbps)
precision = Nb/Nbps;
means = zeros(length(df), precision);
stdv = zeros(length(df), precision);
step = 20;
figure('Name','Gardner')
hold on
ii =1;
Nbexp = 25;
for jj = 1:Nbexp
    % Creation of a random sequence of bits for the moment (row vector)
    bit_tx = randi([0 1],Nb, 1);
    % TX side :
    [symb_tx,signal_tx] = TX(bit_tx, filter,Nbps(j), M);
    signal_tx_M =  kron(signal_tx,copies);

    % Add noise on the signal :
    signal_rx = noise(signal_tx_M, Eb_No_dB,Fsampling,Nb);
    %signal_rx = signal_tx_M;
    %num = size(signal_rx,1);
    % Add CFO
    t = (0:length(signal_rx)-1)./Fsampling;
    signal_cfo =signal_rx .* kron(exp(1j.* (2*pi*df .* t' + phi*ones(length(signal_rx),1))),ones(1,length(Eb_No_dB)));
    %clear signal_rx;

    % RX side :
    symb = zeros(Nb*M/Nbps,length(Eb_No_dB));
    parfor n = 1:length(Eb_No_dB)
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
    means(ii,:) = means(ii,:) + samplingerror(end,:);
    stdv(ii,:) = stdv(ii,:) + samplingerror(end,:).^2;
end
    
    means(ii,:) = means(ii,:)/Nbexp;
    stdv(ii,:) = sqrt(stdv(ii,:)/Nbexp - means(ii,:).^2);
    current = plot(1:step:length(means),means(ii,1:step:end));
    c = get(current, 'color');
    plot(1:step:length(means),means(ii,1:step:end) + stdv(ii,1:step:end),':','Color',c)
    plot(1:step:length(means),means(ii,1:step:end) - stdv(ii,1:step:end),':','Color',c)
    
    
    
    
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