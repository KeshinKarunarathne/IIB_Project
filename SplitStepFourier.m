function [Eout] = SplitStepFourier(Ein, SpS, Rs, CLambda, fibreParams, stepSize)
    
    % Initializing parameters:
    alpha = fibreParams.alpha; % Attenuation constant in dB/km   
    gamma = fibreParams.gamma; % Non-linear fibre constant
    linkLength = fibreParams.linkLength; % Total length of fibre in metres
    D = fibreParams.D; % Fibre disperson constant in ps/nm/km
    
    totalSteps = linkLength/stepSize; % Total number of times the SSF algorithm is applied

    alpha = alpha*log(10)/10; % Converting alpha from dB/km to km-1

    Eout = Ein; % Initialising the output

    for i=1:totalSteps
        timeDomainDistortion = exp(-1i*stepSize*gamma*abs(Eout).^2)*exp(-alpha*stepSize/2).*Eout;
        for j=1:size(Ein,2)
            Eout(:,j) = CDInsertion(timeDomainDistortion(:,j), SpS, Rs, D, CLambda, stepSize*10^3, 1);
        end
    end
end