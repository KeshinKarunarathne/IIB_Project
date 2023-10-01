function g = RRC(span,SpS,rollOff)

    % Initialisations:
    g = zeros(span*SpS+1,1); 
    k = (-span*SpS/2:span*SpS/2)/SpS;
    i1 = find(k==0); 
    i2 = find(abs(4*rollOff*k)-1==0); 
    i3 = 1:length(k);
    i3([i1 i2]) = [];
    k = k(i3);

    % Singularity in k = 0:
    if ~isempty(i1)
        g(i1) = 1-rollOff + 4*rollOff/pi;
    end

    % Singularity in k = 1/(4*rollOff):
    if ~isempty(i2)
        g(i2) = rollOff/sqrt(2)*((1+2/pi)*sin(pi/(4*rollOff))+...
            (1-2/pi)*cos(pi/(4*rollOff)));
    end

    % Calculating the coefficients for k ~= 0 and k ~= 1/(4*rollOff):
    g(i3) = (sin(pi*k*(1-rollOff)) + 4*rollOff*k.*cos(pi*k*(1+rollOff)))...
    ./(pi*k.*(1-(4*rollOff*k).^2));
    
    % Normalizing the amplitude of the filter:
    g = g/max(g);
end
