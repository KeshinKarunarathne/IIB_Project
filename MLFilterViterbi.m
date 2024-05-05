function [wML] = MLFilterViterbi(M, Delta_nu, Rs, OSNRdB, Es, NPol, N)

    % Block length for the Viterbi & Viterbi algorithm, symbol period, and phase noise variance:
    L=2*N+1 ; Ts = 1/Rs ; Sigma_DeltaTheta2 = 2*pi*Delta_nu*Ts;
    
    % Additive noise variance:
    SNRLin=10^(OSNRdB/10)*(2*12.5e9)/(NPol*Rs) ; Sigma_eta2= Es/(2*SNRLin);
    
    % K matrix:
    KAux = zeros(N) ; K = zeros(L);
    for i = 0:N
        for ii = 0:N
            KAux(i+1,ii+1) = min(i,ii);
        end
    end
    K(1:N+1,1:N+1) = rot90(KAux(1:N+1,1:N+1),2);
    K(N+1:L,N+1:L) = KAux(1:N+1,1:N+1);

    % Identity matrix:
    I = eye(L);

    % Obtaining the covariance matrix:
    C = Es^M*M^2*Sigma_DeltaTheta2*K + Es^(M-1)*M^2*Sigma_eta2*I;

    % Filter coefficients:
    wML = (ones(L,1)'/(C)).' ; wML = wML/max(wML);
    
end

