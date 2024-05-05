function [bits, modBits] = SymbolGeneration(modFormat, NSymb)
    % This function maps a set of input bits to the chosen modulation
    % format - QPSK or 16-QAM
    switch modFormat
        case {'QPSK'}
            bitsPerSymbol = 2;
        case '16-QAM'
            bitsPerSymbol = 4;
        otherwise
            error("Modulation format not supported");
    end

    % Generating the seqeunce of bits
    NBits = NSymb*bitsPerSymbol;
    bits = randi([0,1], NBits, 1);
    
    switch modFormat
        case {'QPSK'}
            % In-phase and quadrature bits
            bitsI = bits(2:bitsPerSymbol:NBits); bitsQ = bits(1:bitsPerSymbol:NBits);
    
            % QPSK Modulation {00 -> , 01 -> 10 -> 11 ->}
            xI = 1-2*bitsI; xQ = 1-2*bitsQ;
    
            % Normalised symbols
            modBits = (1/sqrt(2))*(xI + xQ*1i);
        case '16-QAM'
            % In-phase and qaudrature bits
            bitsI1 = bits(4:bitsPerSymbol:NBits); bitsQ1 = bits(3:bitsPerSymbol:NBits);
            bitsI2 = bits(2:bitsPerSymbol:NBits); bitsQ2 = bits(1:bitsPerSymbol:NBits);
    
            % 16-QAM Modulation
            % In-phase
            xI = ((~bitsI2 & ~bitsI1)*(+3) + ( bitsI2 & ~bitsI1)*(+1) + ( bitsI2 & bitsI1)*(-1) + (~bitsI2 & bitsI1)*(-3));
            % Quadrature
            xQ = ((~bitsQ2 & ~bitsQ1)*(+3) + ( bitsQ2 & ~bitsQ1)*(+1) + ( bitsQ2 & bitsQ1)*(-1) + (~bitsQ2 & bitsQ1)*(-3));
    
            % Normalised Symbols
            modBits = (1/sqrt(10))*(xI + xQ*1i);
    end
end