function [E1, E2, E3, E4] = Hybrid(Er, ELo)
    
    % 3-dB coupler transfer function:
    Hc = (1/sqrt(2))*[1 1; 1 -1];

    % ECouplerTL - Signal at the output of the top-left 3-dB coupler at the
    % 90 degree hybrid;
    ECouplerTL = (Hc*[Er.' ; zeros(1,length(Er))]).';

    % ECouplerBL - Signal at the output of the bottom-left 3-dB coupler at
    % the 90 degree hybrid;
    ECouplerBL = (Hc*[ELo.' ; zeros(1,length(Er))]).';

    % ECouplerTR - Signal at the output of the top-right 3-dB coupler at
    % the 90 degree hybrid;
    ECouplerTR = (Hc*[ECouplerTL(:,1).' ; ECouplerBL(:,1).']).';

    % ECouplerBR - Signal at the output of the bottom-right 3-dB coupler at
    % the 90 degree hybrid;
    ECouplerBR = (Hc*[ECouplerTL(:,2).';ECouplerBL(:,2).'*exp(1i*pi/2)]).';
    
    % Output signals:
    E1 = ECouplerTR(:,1); E2 = ECouplerTR(:,2); E3 = ECouplerBR(:,1);
    E4 = ECouplerBR(:,2);
end