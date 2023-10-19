function [E] = LaserModified(Pcw_dBm, centreLambda, lineWidth, SpS, Rs, NSymb, NPol)

    c = 299792458; % Speed of light in m/s

    PcwLinear = 10^(Pcw_dBm/10 - 3); % Linear power of laser output in W

    T = 1/(SpS*Rs); % Period between samples
    t= [0:1:NSymb*SpS - 1]*T; t=t'; % Salient time indices to sample E field
    
    % Generating electric field of chosen wavelength
    if NPol == 1
        E = sqrt(PcwLinear)*exp(1i*(2*pi*c/centreLambda)*t);
    elseif NPol == 2
        E(:,1) = sqrt(PcwLinear)*exp(1i*(2*pi*c/centreLambda)*t);
        E(:,2) = E(:,1);
    else
        error('The possible number of polarizations used must be 1 or 2');
    end
    
    % If the laser linewidth is not 0 Hz, phase noise is inserted 
    if lineWidth ~= 0
        % Calculating the phase noise:
        Var = 2*pi*lineWidth*T ;
        Delta_theta = sqrt(Var)*randn(size(E,1),1);
        Theta = cumsum(Delta_theta);

        % Adding phase noise to the optical signal:
        E = E.*repmat(exp(1i*Theta),1,size(E,2));
    end
end
