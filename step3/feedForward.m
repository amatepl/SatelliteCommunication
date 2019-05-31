function [est_pil, est_cfo] = feedForward(y,pilot,Tsymb,K)
    N = length(pilot);
    L = length(y);
    D = zeros(K, L-N+1);

    for k = 1 : K
        for n = 1 : L - N + 1
            tmp = 0;
            for l = k : N-1
                tmp1 = conj(y(n+l)) * pilot(l+1);
            
                tmp2 = conj(y(n+l-k)) * pilot(l-k+1);
            
                tmp = tmp + tmp1 * conj(tmp2);
            end
            D(k,n) = tmp;
        end 
        D(k,:) = D(k,:)/(N-k);
    end

% Time of arrival estimate
tmp = sum(abs(D),1);
[~, est_pil] = max(tmp);

k = (1:K).' ;
est_cfo = sum(angle(D(k,est_pil))./(2*pi*k*Tsymb), 1);
est_cfo = - est_cfo/K;

end