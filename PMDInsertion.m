function [EOutput, varargout] = PMDInsertion(EInput, DGDSpec, L, N, Rs, SpS, EvalGroupDelay)

    % This function is part of the book Digital Coherent Optical Systems;
    % Darli A. A. Mello and Fabio A. Barbosa;

    % Standard deviation of the Maxwellian distribution:
    SDTau = sqrt(3*pi/8)*DGDSpec;

    % DGD per section (it is equal to the standard deviation per section):
    Tau = (SDTau*sqrt(L*1e-3)/sqrt(N))*1e-12;

    % Frequency vector:
    w = 2*pi*fftshift(-1/2:1/size(EInput,1):1/2-1/size(EInput,1)).'*SpS*Rs;
    
    % Random unitary matrices V and U that describes mode coupling:
    for i = 1:N
        [V(:,:,i),~,U(:,:,i)] = svd(randn(2) + 1i*randn(2));
    end

    % Estimation of the group delay (GD operator):
    if EvalGroupDelay
        % Frequencies to consider when evaluating the group delay:
        wGD = [1 1.1];

        % Obtaining the transfer matrix H:
        for k = 1:numel(wGD)
            % Auxiliary matrix:
            HAux = 1;

            % Delay matrix (Lambda):
            Lambda = [exp(1i*wGD(k)*Tau/2) 0; 0 exp(-1i*wGD(k)*Tau/2)];

            % Transfer function of the i-th section and matrix H:
            for i = 1:N
                Hi = V(:,:,i)*Lambda*U(:,:,i)' ; HAux = HAux*Hi;
            end

            % Matrix H for frequency indicated by k:
            H(:,:,k) = HAux;
        end

        % Obtaining the eigenvalues (num. differentiating H):
        HDiff = (H(:,:,2)-H(:,:,1))/(wGD(2)-wGD(1));
        Eigenvalues = eig(1i*HDiff/(H(:,:,2)));

        % Group delays in ps:
        GroupDelay = abs(real(2*Eigenvalues(1)*1e12));
        varargout{1} = GroupDelay;
    end

    % Signals in frequency domain:
    Freq_E_V = fft(EInput(:,1)) ; 
    Freq_E_H = fft(EInput(:,2));

    for i = 1:N
        % Hermitian of matrix U:
        UHermitian = U(:,:,i)';

        % Rotating the signals according to U:
        E_1 = UHermitian(1,1)*Freq_E_V + UHermitian(1,2)*Freq_E_H;
        E_2 = UHermitian(2,1)*Freq_E_V + UHermitian(2,2)*Freq_E_H;

        % Applying DGD:
        E_1 = exp(1i*w*Tau/2).*E_1 ; E_2 = exp(-1i*w*Tau/2).*E_2;
    
        % Rotating the signals according to V:
        Freq_E_V = V(1,1,i)*E_1 + V(1,2,i)*E_2;
        Freq_E_H = V(2,1,i)*E_1 + V(2,2,i)*E_2;
    end

    % Signals in time domain:
    EOutput(:,1) = ifft(Freq_E_V) ; EOutput(:,2) = ifft(Freq_E_H);
end

