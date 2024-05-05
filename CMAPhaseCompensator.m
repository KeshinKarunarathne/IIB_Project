function [decodedBits, BER] = CMAPhaseCompensator(modFormat, sourceBits, EqSymbs, phaseOffset)
    
    % Initialising parameters of the chosen modulation format
    if strcmp(modFormat, "BPSK")
        M = 2; k = 1;
    elseif strcmp(modFormat, "QPSK")
        M = 4; k = 2;
    else
        error("Entered modulation format not supported. Choose between QPSK and BPSK")
    end

    % Computing initial BER for the two orthoganal polarizations
    decodedInts = pskdemod(EqSymbs, M, phaseOffset);
    decodedBits = int2bit(decodedInts, k);
    BER = ones(2,1);
    for i=1:2
        [~, BER(i)] = biterr(decodedBits(:,i), sourceBits(:,i));
    end
    
    % Finding the phase offset and the correspionding decodedBits which
    % provide a 0 BER for the vertical polarization
    for rotV=1:7
        if (BER(1) > 0.1)
            tempSymbsV = EqSymbs(:,1)*exp(1i*rotV*pi/4);
            tempIntsV = pskdemod(tempSymbsV, M, phaseOffset);
            tempBitsV = int2bit(tempIntsV, k);
            [~, BER(1)] = biterr(tempBitsV, sourceBits(:,1));
            decodedBits(:,1) = tempBitsV;
        else
            break
        end
    end
        
    % Finding the phase offset and the correspionding decodedBits which
    % provide a 0 BER for the horizontal polarization
    for rotH=1:7
        if (BER(2) > 0.1)
            tempSymbsH = EqSymbs(:,2)*exp(1i*rotH*pi/4);
            tempIntsH = pskdemod(tempSymbsH, M, phaseOffset);
            tempBitsH = int2bit(tempIntsH, k);
            [~, BER(2)] = biterr(tempBitsH, sourceBits(:,2));
            decodedBits(:,2) = tempBitsH;
        else
            break
        end
    end

end