function [y,pArr,w11_Arr,w12_Arr,w21_Arr,w22_Arr,errorArr] = AdapEqualizerAlamoutiOriginal(x, trainSymbs, SpS, NTaps, Mu_T, Mu_DD, Mu_p, N1, modFormat)

    % This function implements the 2x2 butterly MIMO for adaptive equalization

    % Inputs
    % x - Input symbols as a single-column vector
    % trainSymbs - original source symbols as a single coulmn vector
    % SpS - samples per symbol
    % NTaps - Number of taps per filter
    % Mu_T - step-size for the training stage
    % Mu_DD - step-size for the decision-directed stage
    % Mu_p - step-size for phase correction
    % N1 - Number of iterations of the training stage before switching to decision-directed mode
    % modFormat - Modulation format -> 'QPSK' for now
    % runs - No. of runs of the equalizer

    % Outputs
    % y - Equalized symbols as a single-column vector
    % pArr - variation of p as a single-column vector
    % w11_Arr - variation of coefficient of w11 
    % w12_Arr - variation of coefficient of w12 
    % w21_Arr - variation of coefficient of w21 
    % w22_Arr - variation of coefficient of w22 
    % errorArr - variation of error as a two-column matrix

    u_o = zeros(length(x)/2,1); u_e = zeros(length(x)/2,1);

    % Separation of odd and even symbols
    count=1;
    for i=1:2*SpS:length(x)-SpS
        u_o(count:count+SpS-1,1) = x(i:i+SpS-1);
        u_e(count:count+SpS-1,1) = x(i+SpS:i+2*SpS-1);
        count = count+SpS;
    end

    % Output length:
    OutLength = size(u_o, 1)/SpS;

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

    % Initializing output symbol streams
    y = zeros(OutLength*2,1);
    
    % Intialising vectors to store data
    pArr = zeros(OutLength,1); vArr = zeros(OutLength,2); errorArr = zeros(OutLength,2);
    w11_Arr = zeros(NTaps, OutLength); w12_Arr = zeros(NTaps, OutLength);
    w21_Arr = zeros(NTaps, OutLength); w22_Arr = zeros(NTaps, OutLength);

    % Initialising filters:
    w11 = zeros(NTaps,1); w12 = zeros(NTaps,1); 
    w21 = zeros(NTaps,1); w22 = zeros(NTaps,1);

    % Initialising filter coefficients
    w11(floor(NTaps/2)+1) = 0.5; w11(floor(NTaps/2)) = 0.5;
    w22(floor(NTaps/2)+1) = 0.5; w22(floor(NTaps/2)) = -0.5; 

    % Initialising phase corrections
    p1a = 1; p1b = 1; p1 = 1; 

    trainSymbIndex = 1;
    for i = 1:OutLength

        % Calculating the outputs:
        v_o = p1*w11.'*u_o(:,i) + conj(p1)*w12.'*conj(u_e(:,i));
        v_e = p1*w21.'*u_o(:,i) + conj(p1)*w22.'*conj(u_e(:,i));

        % Storing equalised symbols for plotting
        vArr(i,1) = v_o; vArr(i,2) = v_e;

        % Error calculation
        if (i <= N1)
            Mu = Mu_T; % Step size for training mode
            e_o = trainSymbs(trainSymbIndex) - v_o;
            e_e = trainSymbs(trainSymbIndex+1) - v_e;
            trainSymbIndex = trainSymbIndex + 2;
        else
            Mu = Mu_DD; % Step size for decision-directed mode
            e_o = DecisionMod(v_o, modFormat) - v_o;
            e_e = DecisionMod(v_e, modFormat) - v_e;
        end

        % Storing errors for plotting
        errorArr(i,1) = e_o; errorArr(i,2) = e_e;

        p1a = p1a + Mu_p*e_o*conj(w11.'*u_o(:,i));
        p1b = p1b + Mu_p*e_o*conj(w12.'*conj(u_e(:,i)));
        p1 = (p1a + conj(p1b))/2; % Phase correction
        pArr(i,1) = p1; % Storing phase correction for plotting

        % Updating the filter coefficients:
        w11 = w11 + Mu*abs(p1)*e_o*conj(u_o(:,i))/p1;
        w12 = w12 + Mu*abs(p1)*e_o*u_e(:,i)/conj(p1);
        w21 = w21 + Mu*abs(p1)*e_e*conj(u_o(:,i))/p1;
        w22 = w22 + Mu*abs(p1)*e_e*u_e(:,i)/conj(p1);
        
        % Storing updated filter coefficients for plotting
        w11_Arr(:,i)=w11; w12_Arr(:,i)=w12; w21_Arr(:,i)=w21; w22_Arr(:,i)=w22;
    
    end

    % Interleaving odd and even symbols
    count = 1;
    for i=1:2:2*OutLength
        y(i,1) = vArr(count,1);
        y(i+1,1) = vArr(count,2);
        count = count + 1;
    end

end