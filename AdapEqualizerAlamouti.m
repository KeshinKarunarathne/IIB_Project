function [y,pArr,w11_Arr,w12_Arr,w21_Arr,w22_Arr] = AdapEqualizerAlamouti(x, trainSymbs, SpS, NTaps, Mu, Mu_p, N1)

    % Serial -> Parallel conversion
    u_o = x(1:2:end, 1); % Odd symbols
    u_e = x(2:2:end, 1); % Even symbols

    % ---------------------------------------------------------------------

    % Rearranging input blocks:
    % Odd symbols
    u_o = [u_o(end-floor(NTaps/2)+1:end,:); u_o ; u_o(1:floor(NTaps/2),:)];
    u_o = convmtx(u_o.',NTaps);
    u_o = u_o(:,NTaps:SpS:end-NTaps+1);

    % Even symbols
    u_e = [u_e(end-floor(NTaps/2)+1:end,:); u_e ; u_e(1:floor(NTaps/2),:)];
    u_e = convmtx(u_e.',NTaps);
    u_e = u_e(:,NTaps:SpS:end-NTaps+1);

    % ---------------------------------------------------------------------

    % Output length:
    OutLength = floor((size(u_o,1)-NTaps+1)/2);

    % Initializing output symbol streams
    v_o = zeros(OutLength,1) ; v_e = zeros(OutLength,1); y = zeros(OutLength*2,1);
    pArr = zeros(OutLength,1);

    % Initial filter coefficients:
    w11 = zeros(NTaps,1); w12 = zeros(NTaps,1); 
    w21 = zeros(NTaps,1); w22 = zeros(NTaps,1);

    p1 = 0; p2 = 0; % Initial phase corrections

    for i = 1:OutLength

        p = (p1 + conj(p2))/2; % Phase correction
        pArr(i) = p;

        % Calculating the outputs:
        v_o(i) = p*w11'*u_o(:,i) + conj(p)*w12'*conj(u_e(:,i));
        v_e(i) = p*w21'*u_o(:,i) + conj(p)*w22'*conj(u_e(:,i));
    
        % Error calculation
        if (i <= N1)
            e_o = trainSymbs(i) - v_o(i);
            e_e = trainSymbs(i) - v_e(i);
        else
            e_o = Decision(v_o(i), modFormat) - v_o(i);
            e_e = Decision(v_e(i), modFormat) - v_e(i);
        end

        % Updating phase correction
        p1 = p1 + Mu_p*e_o*conj(w11'*u_o(:,i));
        p2 = p2 + Mu_p*e_o*conj(w12'*conj(u_e(:,i)));

        % Updating the filter coefficients:
        % Constant modulus algorithm:
        [w11,w12,w21,w22] = AlamoutiEqUpdate(w11, w12, w21, w22, p, e_o, e_e, Mu, u_e, u_o);

        w11_Arr(:,i)=w11; w12_Arr(:,i)=w12; w21_Arr(:,i)=w21; w22_Arr(:,i)=w22;
    
    end

    % Interleaving odd and even symbols
    count = 1;
    for i=1:2:2*OutLength
        y(i,1) = v_o(count);
        y(i+1,1) = v_e(count);
        count = count + 1;
    end


end