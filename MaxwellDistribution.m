function [pdf] = MaxwellDistribution(mu, tau_step, tau_max)
    
    tau = [0:tau_step:tau_max];
    b = mu/(2*sqrt(2/pi));
    pdf = sqrt(2/pi)*(1/b^3)*(tau.^2).*exp(-(tau.^2)/(2*(b^2)));

end
