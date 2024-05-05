function [modulatedBits] = PulseShaping(x, SpS, span, rollOff)
    % This function performs pulse shaping using the RRC filer

    % Obtaining the filter transfer function:
    g = RRC(span,SpS,rollOff);

    % Upsampling the symbols for pulse shaping:
    xUpsamp = upsample(x(:,1),SpS);

    % Filtering the upsampled symbols:
    modulatedBits(:,1) = conv(xUpsamp,g,'same');

    % Normalizing the signal to unitary power:
    modulatedBits(:,1) = modulatedBits(:,1)/sqrt(mean(abs(modulatedBits(:,1)).^2));
end


