function [Ad,Bd] = c2d_zoh_expm(Ac,Bc,Ts)
% ZOH discretisation via Van Loan's expm method (doi:10.1109/tac.1978.1101743).
% Equivalent to c2d(ss(Ac,Bc,eye(n),0), Ts) but avoids the ss/c2d overhead.

    n  = size(Ac,1);
    m  = size(Bc,2);
    M  = expm([Ac Bc; zeros(m,n+m)] * Ts);
    Ad = M(1:n, 1:n);
    Bd = M(1:n, n+1:end);

end