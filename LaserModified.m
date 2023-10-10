function [Ereal, Ecomplex] = LaserModified(Pcw_dBm, centreLambda, lineWidth, SpS, Rs, NSymb, NPol)

    c = 299792458; % Speed of light in m/s
    PcwLinear = 10^(Pcw_dBm/10); % Linear power of laser output
    T = 1/(SpS*Rs); % Period between samples
    t= [0:1:NSymb*SpS - 1]*T; t=t'; % Salient time indices to sample E field
    
    % If the laser linewidth is not 0 Hz, phase noise is inserted 
    Theta = zeros(size(t));
    if lineWidth ~= 0
        % Calculating the phase noise:
        Var = 2*pi*lineWidth*T ;
        Delta_theta = sqrt(Var)*randn(size(Theta,1),1);
        Theta = cumsum(Delta_theta);
    end
    
    % Generating the electric field of the optical signal:
    if NPol == 1
        Ereal = sqrt(PcwLinear)*sin(Theta + (2*pi*c/centreLambda)*t);
        Ecomplex = sqrt(PcwLinear)*exp(1i*(2*pi*c/centreLambda)*t);
    elseif NPol == 2
        Ereal(:,1) = sqrt(PcwLinear)*sin(Theta + (2*pi*c/centreLambda)*t);
        Ereal(:,2) = Ereal(:,1);
        Ecomplex(:,1) = sqrt(PcwLinear)*exp(1i*(2*pi*c/centreLambda)*t);
        Ecomplex(:,2) = Ecomplex(:,1);
    else
        error('The possible number of polarizations used must be 1 or 2');
    end
end
