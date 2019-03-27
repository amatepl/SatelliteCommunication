Fsymbol = 1e6;
beta = 0.3;
RRCTaps = 101;
Fsampling = 8e8;

Delta_t = 1/Fsampling;
t = (-(RRCTaps -1)/2:(RRCTaps -1)/2)*Delta_t;

stepOffset = (1/RRCTaps)*Fsymbol;
highestFreq = stepOffset*(RRCTaps -1)/2;
freqGrid = linspace(-highestFreq,highestFreq,RRCTaps);


interval_3_1 = freqGrid(abs(freqGrid)> (1+beta)*Fsymbol/2  & freqGrid<0);
H1 = zeros(size(interval_3_1));
H = H1;
interval_2_1 = freqGrid((1-beta)*Fsymbol/2 < abs(freqGrid) & abs(freqGrid)<= (1+beta)*Fsymbol/2 & freqGrid<0);
H1 = (1 + cos(pi*(abs(interval_2_1)-(1-beta)*Fsymbol/2)/(Fsymbol*beta)))/(Fsymbol*2);
H = [H,H1];
interval_1_1 = freqGrid(abs(freqGrid)<= (1-beta)*Fsymbol/2 & freqGrid<0);
H1 = ones(size(interval_1_1))/Fsymbol;
H = [H,H1];

interval_1_2 = freqGrid(abs(freqGrid)<= (1-beta)*Fsymbol/2 & freqGrid>=0);
H1 = ones(size(interval_1_2))/Fsymbol;
H = [H,H1];
interval_2_2 = freqGrid((1-beta)*Fsymbol/2 < abs(freqGrid) & abs(freqGrid)<= (1+beta)*Fsymbol/2 & freqGrid>=0);
H1 = (1 + cos(pi*(abs(interval_2_2)-(1-beta)*Fsymbol/2)/(Fsymbol*beta)))/(Fsymbol*2);
H = [H,H1];
interval_3_2 = freqGrid(abs(freqGrid)> (1+beta)*Fsymbol/2  & freqGrid>=0);
H1 = zeros(size(interval_3_2));
H = [H,H1];

h = ifft(H);
h = ifftshift(h);   
size(h)
h = abs(h)/max(abs(h));
figure('Name','h(t)');
plot(t,h)
figure('Name','H(f)');
plot(freqGrid,H)

