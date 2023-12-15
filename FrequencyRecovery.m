function [y, Delta_f] = FrequencyRecovery(x, Rs)

    T = 1/Rs; % Sampling period when SpS=1
    
    % Creating frequency vector for FFT calculation
    f = (-1/2+1/length(x):1/length(x):1/2)*Rs ;

    % Obtaining the FFT for input signal, x, raised to the power of 4
    x_FFT4 = fftshift(abs(fft(x(:,1).^4)));

    % Evaluating frequency offset from FFT above
    Delta_f = (1/4)*f(x_FFT4 == max(x_FFT4));

    % Frequency offset compensation
    k = repmat((0:length(x)-1).',1,size(x,2)); % Temporal indices of x
    y = x.*exp(-1i*2*pi*Delta_f*T*k); % Offset compensation

end