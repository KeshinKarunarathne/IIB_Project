function [v,varargout] = ViterbiCPR(z, Delta_nu, Rs, OSNRdB, Es, NPol, M, ParamViterbi)

    % Block length for the phase estimation:
    N = ParamViterbi.N ; L = 2*N+1;

    % ML filter:
    wML = MLFilterViterbi(M, Delta_nu, Rs, OSNRdB, Es, NPol,  N);

    % Initializing the vector of phase estimates:
    ThetaML = zeros(size(z,1),NPol);

    % Phase noise estimation:
    for Pol = 1:NPol
        % Input blocks:
        zBlocks = [zeros(floor(L/2),1) ; z(:,Pol) ; zeros(floor(L/2),1)];
        zBlocks=convmtx(zBlocks.',L); zBlocks=flipud(zBlocks(:,L:end-L+1));

        % Note that each column of zBlocks is used independently to
        % generate phase estimates:
    
        ThetaML(:,Pol) = (1/M)*angle(wML.'*(zBlocks.^M))-pi/M;
    end
    clearvars zBlocks;

    % Vector of phase estimates after phase unwrapping:
    ThetaPU = zeros(size(ThetaML,1),NPol);

    % Initial 'previous phase' for unwrapping operation:
    ThetaPrev = zeros(1,NPol);

    % Phase unwrapping:
    for i = 1:size(ThetaML,1)
        % Phase unwrapper:
        n = floor(1/2 + (ThetaPrev - ThetaML(i,:))/(2*pi/M));
        ThetaPU(i,:) = ThetaML(i,:) + n*(2*pi/M); ThetaPrev = ThetaPU(i,:);
    end

    % Phase noise compensation:
    v = z.*exp(-1i*ThetaPU);

    % Estimated phase noise as an output of the function:
    if ParamViterbi.PEstimate
        varargout{1} = ThetaPU;
    end
end
