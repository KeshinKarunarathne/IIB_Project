function [out] = CDInsertion(in, SpS, Rs, D, CLambda, L, NPol)

    % This function is part of the book Digital Coherent Optical Systems;
    % Darli A. A. Mello and Fabio A. Barbosa;

    % Constants:
    c = 3e8; % Speed of light

    % Dispersion:
    D = D*1e-6; % In S.I. units

    % Frequency vector:
    w = 2*pi*(-1/2:1/size(in,1):1/2-1/size(in,1)).'*SpS*Rs;

    % Calculating the CD frequency response:
    G = exp(1i*((D*CLambda^2)/(4*pi*c))*L*(w.^2));

    % Inserting CD to the transmitted signal:
    out(:,1) = ifft(ifftshift(G.*fftshift(fft(in(:,1)))));

    % In the case of pol-mux:
    if NPol == 2
        % Inserting CD to the transmitted signal:
        out(:,2) = ifft(ifftshift(G.*fftshift(fft(in(:,2)))));
    end
end

