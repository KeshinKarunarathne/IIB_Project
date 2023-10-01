function [out] = NoiseInsertion(x, modBits, SNRb_dB, SpS)

    % SNR per bit in linear scale:
    SNRb_Lin = 10^(SNRb_dB/10);

    % AWGN standard deviation:
    StdDev = sqrt(mean(abs(x).^2)*SpS/(2*modBits*SNRb_Lin));

    % AWGN generation:
    n = StdDev*randn(length(x),1) + 1i*StdDev*randn(length(x),1);

    % Inserting noise to the signal:
    out = x + n;
end