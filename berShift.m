clear; clc; close all;
% Simulation parameters
Fsymbol = 2e6;
M = 2*100;
Fsampling = M*Fsymbol;
RRCTaps = 155;
beta = 0.3;
Nbps = 4;
Nb = 60000;

bit_tx = randi([0 1],Nb, 1);

Eb_No_dB = 1:1:16;
BER = zeros(length(Nbps),length(Eb_No_dB));

shift=20;   % 0, 4, 10, 20 (0, 0.02, 0.05, 0.1)
% Nyquist
filter = RRCFilter(beta,Fsymbol,Fsampling, RRCTaps);
figure;

    % TX side :
    [symb_tx,signal_tx] = TX(bit_tx, filter,Nbps,M);
    for i = 1:length(Eb_No_dB)
        % Add noise on the signal :
        signal_rx = noise(signal_tx, Eb_No_dB(i),Fsampling,Nb);
        
        % RX side :
        up_symb_rx = conv(signal_rx,filter','full');
        symb = up_symb_rx(RRCTaps:end-RRCTaps+1);
    
        
         %symbrx=zeros(length(symb),1);
         %symbrx(1+shift:end)=symb(1+shift:end,1);

        symb_RX =  downsample(symb(1+shift:end,1),M);
        %symb_RX=downsample(symb,M);
        
        if Nbps > 1 
            bit_rx = demapping(symb_RX, Nbps, 'qam');
        else 
            bit_rx = demapping(symb_RX, Nbps, 'pam');
        end  

        errors = abs(bit_rx - bit_tx);
        BER(i) = sum(errors)/Nb;
    end
    semilogy(Eb_No_dB,BER);
    hold on;

title('BER');

xlabel('Eb/N0 (dB)');
ylabel('BER sampling shift');
grid on;

berTheo2 = berawgn(Eb_No_dB,'qam',2^Nbps);
semilogy(Eb_No_dB,berTheo2,'o');
hold on;
ylim([10e-5 1]);


    