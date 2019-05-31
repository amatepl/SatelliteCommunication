function [corrected, epsilon] = gardner(samples,k,ratio)

% ratio:     Fsample/Fsymbol
% samples:   symb_rx(ii,:) where ii = 1:num, num = size(message_noisy,1) 
% epsilon:   time error ML estimate

% ratio = 1;

size_samples = size(samples);
epsilon = zeros(size_samples(2),ceil(size_samples(1)/ratio));
corrected = zeros(size(epsilon));
samples = [samples; zeros(1,size_samples(2))];          % 0 is added at the and to allow extrapolation.

parfor j = 1:size_samples(2)
    l_epsilon = zeros(1,ceil(size_samples(1)/ratio));
    l_corrected = zeros(1,ceil(size_samples(1)/ratio));
    l_corrected(1) = samples(1,j);
    for n=1:length(l_epsilon)-1

        % Signal linear (here cubic is used) interpolation between two samples.
        % interpolate  =    [y(n - 0.5) y(n)]
        % corrected(n) =    y(n-1)

        interpolate = interp1(1:ratio+1,samples(ratio*(n-1)+1:ratio*n+1,j),[ratio/2+1 ratio+1]- l_epsilon(n),'linear','extrap');
        l_corrected(n+1) = interpolate(2);
        l_epsilon(n+1) = l_epsilon(n) + 2*k*real(interpolate(1)*(conj(l_corrected(n+1)) - conj(l_corrected(n))));
        
    end
    corrected(j,:) = l_corrected;
    epsilon(j,:) = l_epsilon;
end
epsilon = epsilon./ratio;







