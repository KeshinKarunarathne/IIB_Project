function [decodedSequence] = Decision(r, modFormat, bitsOutput)
    
    % Decision:
    switch modFormat
        case 'QPSK'
        % Decision regions for the in-phase components:
        R1 = real(r) >= 0; R2 = real(r) < 0;
        % Applying the decision regions to the imaginary axis:
        R3 = imag(r) >= 0; R4 = imag(r) < 0;
    case '16QAM'
        % Applying the decision regions to the real axis:
        R1 = real(r) >= 2/sqrt(10); R2 = real(r) >= 0;
        R3 = real(r) < 0 ; R4 = real(r) <=-2/sqrt(10);
        % Applying the decision regions to the imaginary axis:
        R5 = imag(r) >= 2/sqrt(10); R6 = imag(r) >= 0;
        R7 = imag(r) < 0 ; R8 = imag(r) <=-2/sqrt(10);
    otherwise
        error('The Supported modulation formats are QPSK and 16-QAM;');
    end

    if bitsOutput
        % Binary labeling:
        switch modFormat
            case 'QPSK'
                ModBits = 2;
                decodedSequence = NaN(length(r),2);
                % Assigning the bits based on the mapping done in the
                % transmitter:
                decodedSequence(R1,2) = zeros(1,sum(R1));
                decodedSequence(R2,2) = ones(1,sum(R2));
                decodedSequence(R3,1) = zeros(1,sum(R3));
                decodedSequence(R4,1) = ones(1,sum(R4));
            case '16QAM'
                ModBits = 4;
                decodedSequence = NaN(length(r),4);
                % Assigning the bits based on the mapping done in the
                % transmitter:
                decodedSequence(R1,[2 4]) = repmat([0 0],sum(R1),1);
                decodedSequence(R2&~R1,[2 4]) = repmat([1 0],sum(R2&~R1),1);
                decodedSequence(R3&~R4,[2 4]) = repmat([1 1],sum(R3&~R4),1);
                decodedSequence(R4,[2 4]) = repmat([0 1],sum(R4),1);
                decodedSequence(R5,[1 3]) = repmat([0 0],sum(R5),1);
                decodedSequence(R6&~R5,[1 3]) = repmat([1 0],sum(R6&~R5),1);
                decodedSequence(R7&~R8,[1 3]) = repmat([1 1],sum(R7&~R8),1);
                decodedSequence(R8,[1 3]) = repmat([0 1],sum(R8),1);
        end

        % Obtaining the decided bits as a column vector:
        decodedSequence = reshape(decodedSequence',1,length(r)*ModBits)';

    else
        % decodedSequence Symbols:
        decodedSequence = zeros(size(r));
        switch modFormat
            case 'QPSK'
                % Assigning the bits based on the mapping done in the
                % transmitter:
                decodedSequence(R1) = decodedSequence(R1) + (1/sqrt(2));
                decodedSequence(R2) = decodedSequence(R2) - (1/sqrt(2));
                decodedSequence(R3) = decodedSequence(R3) + (1i*1/sqrt(2));
                decodedSequence(R4) = decodedSequence(R4) - (1i*1/sqrt(2));
            case '16QAM'
                % Assigning the bits based on the mapping done in the
                % transmitter:
                decodedSequence(R1) = decodedSequence(R1) + (3/sqrt(10));
                decodedSequence(R2 & ~R1) = decodedSequence(R2&~R1) + (1/sqrt(10));
                decodedSequence(R3 & ~R4) = decodedSequence(R3&~R4) - (1/sqrt(10));
                decodedSequence(R4) = decodedSequence(R4) - (3/sqrt(10));
                decodedSequence(R5) = decodedSequence(R5) + (3i/sqrt(10));
                decodedSequence(R6 & ~R5) = decodedSequence(R6&~R5) + (1i/sqrt(10));
                decodedSequence(R7 & ~R8) = decodedSequence(R7&~R8) - (1i/sqrt(10));
                decodedSequence(R8) = decodedSequence(R8) - (3i/sqrt(10));
        end
    end
end
