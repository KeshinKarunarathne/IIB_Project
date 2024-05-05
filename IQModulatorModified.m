function [Eout] = IQModulatorModified(vI, vQ, Ein, V_pi)

    phiI = vI*pi/V_pi; phiQ = vQ*pi/V_pi;

    Eout = (0.5*cos(0.5*phiI) + 0.5i*cos(0.5*phiQ)).*Ein;

    % Debugging line modelling a single MZ modulator
    % Eout = Ein.*cos(phiI/2);
end