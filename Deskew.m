function [rOut] = Deskew(rIn, SpSRx, Rs, NPol, N, ParamSkew)

    % Skew to be compensated using the TADC as reference:
    TADC=1/(SpSRx*Rs) ; 
    Skew = [ParamSkew.TauIV/TADC ParamSkew.TauQV/TADC];

    if NPol == 2
        Skew = [Skew ParamSkew.TauIH/TADC ParamSkew.TauQH/TADC];
    end

    % Using the min skew as reference:
    Skew = Skew - min(Skew);

    % Integer and fractional part of the skew:
    nTADC = floor(Skew) ; muTADC = -(Skew-nTADC);

    % Obtaining the FIR filter and interpolating the signals:
    NTaps = N+1; % Number of filter taps;
    for i = 1:size(rIn,2)
        L = zeros(NTaps,1) ; Aux = 1;

        % Obtaining the Lagrangean interpolator:
        for n = (0:N) - floor(mean(0:N)) + nTADC(i)
            m = (0:N) - floor(mean(0:N)) + nTADC(i) ; m(m == n) = [];
            L(Aux) = prod((muTADC(i) - m)./(n - m)) ; Aux = Aux + 1;
        end

        % Interpolating the received signal (sIn):
        sAux = flipud(convmtx([zeros(1,floor(NTaps/2)) rIn(:,i).'...
            zeros(1,floor(NTaps/2))],NTaps));
        sAux = sAux(:,NTaps:end-(NTaps)+1) ; 
        rOut(:,i) = (L.'*sAux).';
    end
end