function [v, varargout] = BPS(z, modOrder, NPol, ParamBPS)

    % This function is part of the book Digital Coherent Optical Systems;
    % Darli A. A. Mello and Fabio A. Barbosa;

    % Number of test rotations and ambiguity angle; Total BPS block length:
    B = ParamBPS.B; p = pi/2 ; N = ParamBPS.N; L=2*N+1;

    % Creating a vector of test carrier phase angles:
    b = -B/2:1:B/2-1 ; ThetaTest = p*b/B;

    % Creating a matrix of test carrier phase angles:
    ThetaTestMatrix = repmat(exp(-1j*ThetaTest),L,1);
    if NPol == 2
        ThetaTestMatrix = cat(3,ThetaTestMatrix,ThetaTestMatrix);
    end

    % Input blocks:
    % V-pol. orientation:
    zB_V = [zeros(floor(L/2),1); z(:,1) ; zeros(floor(L/2),1)];
    zB_V = convmtx(zB_V.',L) ; zB_V = flipud(zB_V(:,L:end-L+1));

    if NPol == 2
        % H-pol. orientation:
        zB_H = [zeros(floor(L/2),1) ; z(:,2) ; zeros(floor(L/2),1)];
        zB_H = convmtx(zB_H.',L) ; zB_H = flipud(zB_H(:,L:end-L+1));
        zBlocks = cat(3,zB_V,zB_H) ; clearvars zBlocks_V zBlocks_H;
    else
        zBlocks = zB_V; clearvars zBlock_V;
    end

    % Vector of phase estimates and initial phase for the PU:
    ThetaPU = zeros(size(zBlocks,2),NPol) ; ThetaPrev = zeros(1,NPol);

    % Phase noise estimates:
    for i = 1:size(zBlocks,2)
        % Applying the test phase angles to the symbols:
        zRot = repmat(zBlocks(:,i,:),1,B,1).*ThetaTestMatrix;

        % Decision of the rotated symbols:
        zRot_Decided = qammod(qamdemod(zRot, modOrder, UnitAveragePower=true), modOrder, UnitAveragePower=true); % 16-QAM
        % zRot_Decided = Decision(zRot,ModFormat,false);

        % Intermidiate signal to be minimized:
        m = sum(abs(zRot-zRot_Decided).^2,1);

        % Estimating the phase noise as the angle that minimizes 'm':
        [~,im] = min(m,[],2) ; Theta = reshape(ThetaTest(im),1,NPol);

        % Applying the phase unwrapper to the estimated phase angles:
        ThetaPU(i,:) = Theta + floor(0.5 - (Theta-ThetaPrev)./(p)).*(p);

        % Updating the previous phase variable:
        ThetaPrev = ThetaPU(i,:);
    end

    % Compensating for the phase noise:
    v = z.*exp(-1i*ThetaPU);

    % Estimated phase noise as an output of the function:
    if ParamBPS.PEstimate
        varargout{1} = ThetaPU;
    end
end