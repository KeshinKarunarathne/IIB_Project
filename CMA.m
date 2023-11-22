function [w1V,w1H,w2V,w2H] = CMA(xV, xH, y1, y2, w1V, w1H, w2V, w2H, R, Mu)
    
    % Updating the filters:
    w1V = w1V + Mu*xV*(R-abs(y1).^2)*conj(y1);
    w1H = w1H + Mu*xH*(R-abs(y1).^2)*conj(y1);
    w2V = w2V + Mu*xV*(R-abs(y2).^2)*conj(y2);
    w2H = w2H + Mu*xH*(R-abs(y2).^2)*conj(y2);

end
