function [EOutput, varargout] = PMDInsertionMod(EInput, DGDSpec, L, N, Rs, SpS)

    % Standard deviation of the Maxwellian distribution:
    SDTau = sqrt(3*pi/8)*DGDSpec;

    % DGD per section (it is equal to the standard deviation per section):
    Tau = (SDTau*sqrt(L*1e-3)/sqrt(N))*1e-12;

    % Frequency vector:
    % w = 2*pi*fftshift(-1/2:1/size(EInput,1):1/2-1/size(EInput,1)).'*SpS*Rs;
    
    % Random unitary matrices V and U that describes mode coupling:
    % for i = 1:N
        % [V(:,:,i),~,U(:,:,i)] = svd(randn(2) + 1i*randn(2));
    % end

    % Signals in frequency domain:
    % Freq_E_V = fft(EInput(:,1)); 
    % Freq_E_H = fft(EInput(:,2));

    EOutput = EInput;

    for i = 1:N
        theta = randi(4,1)*180/8;

        EOutput = RotationMatrix(EOutput, theta);

        % Applying DGD:
        EOutput(:,1) = exp(1i*w*Tau/2).*E_1 ; E_2 = exp(-1i*w*Tau/2).*E_2;
    
        % Rotating the signals according to V:
        Freq_E_V = V(1,1,i)*E_1 + V(1,2,i)*E_2;
        Freq_E_H = V(2,1,i)*E_1 + V(2,2,i)*E_2;
    end

end

