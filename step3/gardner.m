function [corrected, epsilon] = gardner(samples,k,ratio)

% ratio:     Fgardner/Fsymbol with Fgardner = Fsampling/2
% samples:   symb_rx(ii,:) where ii = 1:num, num = size(message_noisy,1) 
% epsilon:   time error ML estimate

epsilon = zeros(1,ceil(length(samples)/ratio));
corrected = zeros(size(epsilon));
corrected(1) = samples(1);

for n=1:length(epsilon)-1
    
    % Signal linear (here cubic is used) interpolation between two samples.
    % interpolate  =    [y(n - 0.5) y(n)]
    % corrected(n) =    y(n-1)
    
    interpolate = interp1(1:ratio+1,samples(ratio*(n-1)+1:ratio*n+1),[ratio/2+1 ratio+1],'linear');
    
    corrected(n+1) = interpolate(2);
    epsilon(n+1) = epsilon(n) + 2*k*real(interpolate(1)*(conj(corrected(n+1)) - conj(corrected(n))));
    %disp(num2str(epsilon(n+1)));
    %disp(num2str(2*k*real(interpolate(1)*(conj(corrected(n+1)) - conj(corrected(n))))));
end

epsilon = epsilon./ratio;







