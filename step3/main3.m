clear; clc; close all;
addpath('../step1');
addpath('../shared');

% Simulation parameters
Fsymbol = 2e6;
M = 2;
Fsampling = M*Fsymbol;
RRCTaps = 155;
beta = 0.3;
Nbps = [1];
precision = 1e6;
k = 0.01;
Fgardner = Fsampling/2;

fc=2e9;

Nb = 60000;
% Creation of a random sequence of bits for the moment (row vector)
bit_tx = randi([0 1],Nb, 1);

Eb_No_dB = 1:2:16;
copies = ones(1,length(Eb_No_dB));
BER = zeros(length(Nbps),length(Eb_No_dB));

%df = (1e-6*fc)*25; % ppm
df= 0;
phi=pi/16;
% Nyquist
filter = RRCFilter(beta,Fsymbol,Fsampling, RRCTaps);
figure;
j=1;
%for j = 1:length(Nbps)
    % TX side :
    [~,signal_tx] = TX(bit_tx, filter,Nbps(j), M);
    signal_tx_M =  kron(signal_tx,copies);
    %for i = 1:length(Eb_No_dB)
        % Add noise on the signal :
        %   size(signal_rx) = number of bits x length(Eb_No_dB)
        signal_rx = noise(signal_tx_M, Eb_No_dB,Fsampling,Nb);
        %num = size(signal_rx,1);
        num = 1;
        % Add CFO
        t = (0:length(signal_rx)-1)./Fsampling;
        signal_cfo =signal_rx .* kron(exp(1j.* (2*pi*df .* t')),ones(1,length(Eb_No_dB)));
        %clear signal_rx;
        
        % RX side :
        %up_symb_rx = zeros(size(signal_cfo));
        symb = zeros(Nb*M/Nbps,length(Eb_No_dB));
        parfor n = 1:length(Eb_No_dB)
            up_symb_rx = conv(signal_cfo(:,n),filter','full');
            symb(:,n) = up_symb_rx(RRCTaps:end-RRCTaps+1);
        end
        
        %symb = up_symb_rx(RRCTaps:end-RRCTaps+1,:);
        size_symb = size(symb);
%         t = (0:length(symb)-1)./Fsampling;
        t = (0:size_symb(1)-1)./Fsampling;
        %clear signal_cfo;
        
        symb_rx = downsample(symb,M);
        precision = Nb/Nbps(j);
        %corrected = zeros(num,precision);
        %epsilon = zeros(num,precision);
%         for ii=1:num
%             [corrected(ii,:), epsilon(ii,:)] = gardner(symb_rx(ii,:),k,Fgardner/Fsymbol);
%         end
        [corrected, epsilon] = gardner(symb_rx,k,Fgardner/Fsymbol);
        
        symb_rx2 = corrected';
        bit_rx = zeros(Nb,length(Eb_No_dB));
        parfor i = 1:length(Eb_No_dB)
            l_bit_rx = zeros(length(Eb_No_dB),1);
            if Nbps(j) > 1 
                l_bit_rx = demapping(symb_rx2(:,i), Nbps(j), 'qam');
            else 
                l_bit_rx = demapping(symb_rx2(:,i), 1, 'pam');
            end 
            bit_rx(:,i) = l_bit_rx;
        end
        errors = abs(bit_rx - kron(bit_tx,ones(1,length(Eb_No_dB))));
        BER(j,:) = sum(errors)/Nb;
%    end
    semilogy(Eb_No_dB,BER(j,:));
    hold on;
%end
title('Effect of CFO on BER (ISI Only)');
legend('BPSK','QPSK','QAM16','QAM64');
xlabel('Eb/N0 (dB)');
ylabel('BER');
grid on;
hold on;
berTheo1 = berawgn(Eb_No_dB,'pam',2^Nbps(1));
semilogy(Eb_No_dB,berTheo1,'-o');
hold on;
% berTheo2 = berawgn(Eb_No_dB,'qam',2^Nbps(2));
% semilogy(Eb_No_dB,berTheo2,'-o');
% hold on;
% berTheo3 = berawgn(Eb_No_dB,'qam',2^Nbps(3));
% semilogy(Eb_No_dB,berTheo3,'-o');
% hold on;
% berTheo4 = berawgn(Eb_No_dB,'qam',2^Nbps(4));
% semilogy(Eb_No_dB,berTheo4,'-o');
ylim([10e-5 1]);
%legend('BPSK theo','QPSK theo','QAM16 theo','QAM64 theo');

% [symb_tx,signal_tx] = TX(bit_tx, filter,2, M);
% figure;
% plot(symb_tx,'o');
% hold on;
% plot(cfo(signal_rx,df,phi,Fsampling),'x')
% legend('Emitted symbols', 'phase offset');
% title('Effect of phase offset pi/16')
% grid on;


    