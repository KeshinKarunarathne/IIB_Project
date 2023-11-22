function [k, symbols] = PlotConstellation(modFormat, phaseOffset)
    
    switch modFormat
        % Defining parameters values specific to each modulation format
        case "BPSK"
            modOrder = 2;
        case "QPSK"
            modOrder = 4;
        case "16QAM"
            modOrder = 16;
        otherwise
            error('The supported modulation formats are BPSK, QPSK and 16-QAM');
    end

    k= log2(modOrder); % Bits per symbols
    alphabet = 0:modOrder-1; % Alphabet in integers
    
    switch modFormat
        % Defining the constellation for the chosen modulation format
        case "16QAM"
            symbols = qammod(alphabet, modOrder); 
        case "BPSK"
            symbols = pskmod(alphabet, modOrder, phaseOffset); 
        otherwise
            symbols = pskmod(alphabet, modOrder, phaseOffset); 
    end

    % Plotting symbol constellation
    plot(symbols, ".", MarkerSize=20);
    title(strcat(modFormat, " Constellation"));
    xlabel("Real");
    ylabel("Imaginary");
    grid on;
    hold on;
    % phase = [0:pi/20:2*pi];
    % circleX = cos(phase); circleY = sin(phase);
    % plot(circleX, circleY, "--");
    xlim([-4,4]); ylim([-4,4]);
    hold off

end
