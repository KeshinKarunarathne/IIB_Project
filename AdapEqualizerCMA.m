function [y,w1V_Arr,w1H_Arr,w2V_Arr,w2H_Arr] = AdapEqualizerCMA(x, SpS, NTaps, Mu, SingleSpike, N1, NRemove)

    % This function is part of the book Digital Coherent Optical Systems;
    % Darli A. A. Mello and Fabio A. Barbosa;

    % Radius for the CMA (E{|x|^4}/E{|x|^2}):
    R_CMA = 1; 

    % Input blocks:
    x = [x(end-floor(NTaps/2)+1:end,:) ; x ; x(1:floor(NTaps/2),:)];
    xV = convmtx(x(:,1).',NTaps) ; xH = convmtx(x(:,2).',NTaps);
    xV = xV(:,NTaps:SpS:end-NTaps+1) ; xH = xH(:,NTaps:SpS:end-NTaps+1);

    % Output length:
    OutLength = floor((size(x,1)-NTaps+1)/2) ; clearvars x

    % Initializing the outputs
    y1 = zeros(OutLength,1) ; y2 = zeros(OutLength,1);

    % Initial filter coefficients:
    w1V = zeros(NTaps,1); w1H = zeros(NTaps,1); 
    w2V = zeros(NTaps,1); w2H = zeros(NTaps,1);

    % If single spike initialization:
    if SingleSpike
        w1V(floor(NTaps/2)+1) = 1;
    end

    for i = 1:OutLength

        % Calculating the outputs:
        y1(i) = w1V'*xV(:,i) + w1H'*xH(:,i);
        y2(i) = w2V'*xV(:,i) + w2H'*xH(:,i);

        % Updating the filter coefficients:
        % Constant modulus algorithm:
        [w1V,w1H,w2V,w2H] = CMA(xV(:,i),xH(:,i),y1(i),y2(i),w1V,w1H,w2V,w2H,R_CMA,Mu);

        w1V_Arr(:,i)=w1V; w1H_Arr(:,i)=w1H; w2V_Arr(:,i)=w2V; w2H_Arr(:,i)=w2H;
    
        % Reinitialization of the filter coefficients:
        if i == N1 && SingleSpike
            w2H = conj(w1V(end:-1:1,1)) ; w2V = -conj(w1H(end:-1:1,1));
        end
    end

    % Output samples:
    y = [y1 y2] ; y = y(1+NRemove:end,:);

end
