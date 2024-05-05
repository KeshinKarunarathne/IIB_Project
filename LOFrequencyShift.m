function [EOut] = LOFrequencyShift(EIn, Delta_f, T)

    % Phase offset between consecutive samples due to 'Delta_f':
    DeltaTheta_f = -2*pi*Delta_f*T;

    % Inserting the effects of frequency offset into signal 'EIn':
    k = repmat((0:size(EIn,1)-1).',1,size(EIn,2));
    EOut = EIn.*exp(1i*k*DeltaTheta_f);
    
end