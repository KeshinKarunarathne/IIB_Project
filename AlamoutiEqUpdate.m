function [w11, w12, w21, w22] = AlamoutiEqUpdate(w11, w12, w21, w22, p, e_o, e_e, Mu, u_e, u_o)
    % Updating filter cofficients for adaptive equalizer for Alamouti
    % Coding

    w11 = w11 + Mu*abs(p)*e_o*conj(u_o)/p;
    w12 = w12 + Mu*abs(p)*e_o*u_e/conj(p);
    w21 = w21 + Mu*abs(p)*e_e*conj(u_o)/p;
    w22 = w22 + Mu*abs(p)*e_e*u_e/conj(p);

end