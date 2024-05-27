function [E] = Laser(Pcw, laserLinewidth, SpS, Rs, NSymb, NPol)

    % This function is part of the book Digital Coherent Optical Systems;
    % Darli A. A. Mello and Fabio A. Barbosa;

    % Calculating the linear power of the continuous wave:
    PcwLinear = 1e-3*10^(Pcw/10);

    % Generating the electric field of the optical signal:
    if NPol == 1
        E = ones(SpS*NSymb,1)*sqrt(PcwLinear);
    elseif NPol == 2
        E = ones(SpS*NSymb,2)*sqrt(PcwLinear/2);
    else
        error('The possible number of polarizations used must be 1 or 2');
    end

    % If the laser linewidth is not 0 Hz, phase noise is inserted 
    if laserLinewidth ~= 0
        % Period between samples at the (oversampled) transmit. signal:
        T = 1/(SpS*Rs);
        
        % Calculating the phase noise:
        Var = 2*pi*laserLinewidth*T ;
        Delta_theta = sqrt(Var)*randn(size(E,1),1);
        Theta = cumsum(Delta_theta);

        % Adding phase noise to the optical signal:
        E = E.*repmat(exp(1i*Theta),1,size(E,2));
    end
end
