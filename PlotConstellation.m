function [k, symbols] = PlotConstellation(modOrder, phaseOffset)
    
    if modOrder == 2
        modFormat = "BPSK";
    elseif modOrder == 4
        modFormat = "QPSK";
    else
        error('The supported modulation formats are QPSK and 16-QAM');
    end

    k= log2(modOrder); % Bits per symbols
    alphabet = 0:modOrder-1; % Alphabet in integers

    % Symbols of the chosen modulation format
    symbols = pskmod(alphabet, modOrder, phaseOffset); 

    % Plotting symbol constellation
    plot(real(symbols), imag(symbols), ".", MarkerSize=20);
    title(strcat(modFormat, " Constellation"));
    xlabel("Real");
    ylabel("Imaginary");
    grid on;
    hold on;
    phase = [0:pi/20:2*pi];
    circleX = cos(phase); circleY = sin(phase);
    plot(circleX, circleY, "--");
    xlim([-2,2]); ylim([-2,2]);
    hold off

end
