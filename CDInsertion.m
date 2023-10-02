function [out] = CDInsertion(in, SpS, Rs, D,CLambda, L, NPol)

    % Constants:
    c = 299792458; % Speed of light

    % Dispersion:
    D=D*1e-12/(1e-9*1e3); % In S.I. untits

    % Frequency vector:
    w=2*pi*(-1/2:1/size(in,1):1/2-1/size(in,1)).'*SpS*Rs;

    % Calculating the CD frequency response:
    G = exp(1i*((D*CLambda^2)/(4*pi*c))*L*w.^2);

    % Inserting CD to the transmitted signal:
    out(:,1) = ifft(ifftshift(G.*fftshift(fft(in(:,1)))));

    % In the case of pol-mux:
    if NPol == 2
        % Inserting CD to the transmitted signal:
        out(:,2) = ifft(ifftshift(G.*fftshift(fft(in(:,2)))));
    end
end

