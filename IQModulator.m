function [Eout] = IQModulator(xb, Ein, ParamMZM)

    % This function is part of the book Digital Coherent Optical Systems;
    % Darli A. A. Mello and Fabio A. Barbosa;

    % Obtaining the in-phase and quadrature components of the electrical:
    mI = real(xb); mI = mI/max(abs(mI)); % In-phase;
    mQ = imag(xb); mQ = mQ/max(abs(mQ)); % Quadrature;

    % Setting the signal excursion:
    mI = mI*(ParamMZM.MaxExc-ParamMZM.MinExc)/2;
    mQ = mQ*(ParamMZM.MaxExc-ParamMZM.MinExc)/2;

    % Obtaining the signals after considering the bias:
    vI = mI+ParamMZM.Bias; vQ = mQ+ParamMZM.Bias;

    % Phase modulation in the in-phase and quadrature branches;
    PhiI = pi*(vI)/ParamMZM.Vpi; PhiQ = pi*(vQ)/ParamMZM.Vpi;

    % IQM output signal:
    Eout = (0.5*cos(0.5*PhiI) + 0.5i*cos(0.5*PhiQ)).*Ein;
end
