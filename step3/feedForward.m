function [pilPos deltaF] = feedForward(y,pilot,N,K,fsymb)
    nmax = 200;
    D = zeros(K,nmax);
    y_con = conj(y);
    for k = 1:K
        for l = k:N-1
            for n = 1:nmax
                y1 = y(n+l).*pilot(l+1);
                y2 = conj(y(n+l-k).*pilot(l-k+1));
                D(k,n) = D(k,n)+y1.*y2;
            end
        end
        D(k,:) = D(k,:)/(N-k);
    end
    pilPos = max(sum(abs(D)));
    deltaF = -1/K*sum(fsymb*angle(D(:,pilPos))./(2*pi*(1:K).'));
end