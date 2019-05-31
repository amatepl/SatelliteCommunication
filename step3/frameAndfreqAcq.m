clear; clc; close all;
addpath('../step1');
addpath('../shared');
% Simulation parameters
modulation = 'qam';
Fsymbol = 1e6;
M = 8;

Fsymbol = 2*Fsymbol;
Tsymbol = 1/Fsymbol;
Fsamp = 8*Fsymbol;
Tsamp = 1/Fsamp;

RRCTaps = 365;
beta = 0.3;
filter = RRCFilter(beta,Fsymbol,Fsamp, RRCTaps);

Nbps = [1 2 4 6];
j = 2;
Eb_No_dB = [0:1:16];

fc=2e9;
df = 0;
phi= 0;
shift=0;
pilot_pos = 10;

pilotLength = [10 20 40];
K = [1 8 16];

Nb = 200;

expNum = 1200;

for z = 1:expNum
  for k = 1:length(K)
    for n = 1:length(pilotLength)
      % Creation of a random sequence of bits for the moment (row vector)
      bit_data = randi(2,1,Nb)-1;
      bit_pilot = randi(2,1,pilotLength(n))-1;

      symb_data = mapping(bit_data', Nbps(j), modulation);
      symb_pilot = mapping(bit_pilot', Nbps(j), modulation);
      symb_tx = [symb_data(1:pilot_pos-1); symb_pilot;symb_data(pilot_pos:end)];
      up_symb_tx = upsample(symb_tx,M);
      signal_tx = conv(up_symb_tx,filter');
      
      for i = 1:length(Eb_No_dB)
        % Add noise on the signal :
        signal_rx = noise(signal_tx, Eb_No_dB(i),Fsamp,Nb);
        % Add CFO
        arg1 = ((0:length(signal_tx)-1)-(RRCTaps-1)/2)*Tsamp+phi;
        signal_cfo =signal_rx .* exp(1j.* (2*pi*df*arg1'));

        % RX side :
        up_symb_rx = conv(signal_cfo,filter');
        up_symb_rx = up_symb_rx(RRCTaps:end-RRCTaps+1);
        symb_rx = downsample(up_symb_rx(1+shift:end),M);
        
        arg2 = (0:length(symb_rx)-1)*Tsamp*M;
        symb_cfo = symb_rx.*exp(-1j*(2*pi*df*arg2'));
        [est_pil, est_cfo] = feedForward(symb_rx,symb_pilot,Tsymbol,K(k));
        
        bit_rx = demapping(symb_cfo, Nbps(j), modulation);
        pilot_error(z,i,k,n) = est_pil - pilot_pos;
        cfo_error(z,i,k,n) = (df + est_cfo)*1e6/fc;
      end
    end
  end
end
pilot_mean = std(pilot_error);
cfo_mean = std(cfo_error);

figure();
plot(Eb_No_dB,pilot_mean(1,:,2,1),'-o');
hold on;
plot(Eb_No_dB,pilot_mean(1,:,2,2),'-o');
hold on;
plot(Eb_No_dB,pilot_mean(1,:,2,3),'-o');
title('Pilot ToA error as a function of N for PAM modulation')
xlabel('Eb/No [dB]');
ylabel('Time error stdev [samples]');
legend('N = 10, K = 8','N = 20, K = 8', 'N = 40, K = 8');
grid on;

figure();
plot(Eb_No_dB,pilot_mean(1,:,1,2),'-o');
hold on;
plot(Eb_No_dB,pilot_mean(1,:,2,2),'-o');
hold on;
plot(Eb_No_dB,pilot_mean(1,:,3,2),'-o');
title('Pilot ToA error as a function of K for PAM modulation')
xlabel('Eb/No [dB]');
ylabel('Time error stdev [samples]');
legend('N = 20, K = 1','N = 20, K = 8', 'N = 20, K = 16');
grid on;

figure();
plot(Eb_No_dB,cfo_mean(1,:,2,1),'-o');
hold on;
plot(Eb_No_dB,cfo_mean(1,:,2,2),'-o');
hold on;
plot(Eb_No_dB,cfo_mean(1,:,2,3),'-o');
title('CFO error as a function of N for PAM modulation')
xlabel('Eb/No [dB]');
ylabel('Frequency error stdev[ppm]');
legend('N = 10, K = 8','N = 20, K = 8', 'N = 40, K = 8');
grid on;

figure();
plot(Eb_No_dB,cfo_mean(1,:,1,2),'-o');
hold on;
plot(Eb_No_dB,cfo_mean(1,:,2,2),'-o');
hold on;
plot(Eb_No_dB,cfo_mean(1,:,3,2),'-o');
title('CFO error as a function of K for PAM modulation')
xlabel('Eb/No [dB]');
ylabel('Frequency error stdev[ppm]');
legend('N = 20, K = 1','N = 20, K = 8', 'N = 20, K = 16');
grid on;