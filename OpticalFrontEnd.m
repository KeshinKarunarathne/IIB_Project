function [Out] = OpticalFrontEnd(Er, ELo, R, NPol)
    
    % 90 degree Hybrid:
    [E1,E2,E3,E4] = Hybrid90(Er(:,1),ELo(:,1));

    % Photodetection:
    i1 = R*(E1.*conj(E1)) ; i2 = R*(E2.*conj(E2)); % In-phase;
    i3 = R*(E3.*conj(E3)) ; i4 = R*(E4.*conj(E4)); % Quadrature;

    % Balanced photodetection:
    iI = i1 - i2 ; iQ = i3 - i4;

    % Output signal:
    Out(:,1) = iI ; Out(:,2) = iQ;

    % In the case of pol-mux:
    if NPol == 2
        % 90 degree Hybrid:
        [E1,E2,E3,E4] = Hybrid90(Er(:,2),ELo(:,2));

        % Photodetection:
        i1 = R*(E1.*conj(E1)) ; i2 = R*(E2.*conj(E2)); % In-phase;
        i3 = R*(E3.*conj(E3)) ; i4 = R*(E4.*conj(E4)); % Quadrature;

        % Balanced photodetection:
        iI = i1 - i2 ; iQ = i3 - i4;
        
        % Output signal:
        Out(:,3) = iI ; Out(:,4) = iQ;
    end
end