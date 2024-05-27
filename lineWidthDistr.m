function [lineWidth_Arr, BERQAM_shaped_final] = lineWidthDistr(distr, nuTs_Arr)
    NSymb = 50000; % No. of symbols 
    SpS = 2; % Samples per symbol
    Rs = 50e9; % Symbol rate in symbols/sec
    T = 1/(SpS*Rs); % Sampling period in s
    t= [0:1:NSymb*SpS - 1]*T; % Time vector for plotting time-series data

    %% 
    
    % Generating random bit sequence(s)
    sourceBitsQAM = randi([0 1], NSymb*4, 1); % Matrix of source bits for QAM
    sourceIntsQAM = bit2int(sourceBitsQAM, 4); % Vector of source bits converted to ints
    sourceSymbolsQAM = qammod(sourceIntsQAM, 16, UnitAveragePower=true); % Maps source bit vector to symbol vector
    codedSymbolsQAM = AlamoutiCoding(sourceSymbolsQAM); % Encodes symbols using Alamouti coding
    
    vI_QAM = real(codedSymbolsQAM);
    vQ_QAM = imag(codedSymbolsQAM);
    
    vI_QAM_upsampled = repelem(vI_QAM, SpS, 1); % Upsampling according to SpS
    vQ_QAM_upsampled = repelem(vQ_QAM, SpS, 1); % Upsampling according to SpS

    %% 


    % RRC filter parameters
    filterSpan = 16; % Span of RRC filter in no. of symbols
    rollOff = 0.1;
    filterResponse = rcosdesign(rollOff, filterSpan, SpS, "sqrt")';
    samples = 1+(NSymb-1)*SpS + filterSpan*SpS;

    % ---------------------------------------------------------------------------

    % Using RRC filter to shape vI and vQ
    vI_QAM_shaped = upfirdn(vI_QAM, filterResponse, SpS, 1); vQ_QAM_shaped = upfirdn(vQ_QAM, filterResponse, SpS, 1);

    %% 

    % Initialising laser parameters
    lineWidth_Arr = nuTs_Arr*Rs;
    laserPw_dBm = 30;
    laserPw_Lin = 1e-3*10^(laserPw_dBm/10);

    modOpticSigQAM = zeros(size(vI_QAM_upsampled,1), 2, length(lineWidth_Arr)); 
    modOpticSigQAM_shaped = zeros(size(vI_QAM_shaped,1), 2, length(lineWidth_Arr));
    % Modulating symbols using different laser linewidth
    for i=1:length(lineWidth_Arr)
        % Generating output from laser modelled as a continuous-wave source
        laserE = Laser(laserPw_dBm, lineWidth_Arr(i)*distr, SpS, Rs, NSymb, 2);
        laserE_shaped = Laser(laserPw_dBm, lineWidth_Arr(i)*distr, SpS, Rs, samples/SpS, 2);
    
        modOpticSigQAM(:,:,i) = (vI_QAM_upsampled+1i*vQ_QAM_upsampled).*laserE;
        modOpticSigQAM_shaped(:,:,i) = (vI_QAM_shaped+1i*vQ_QAM_shaped).*laserE_shaped;
    end

    %% 
    
    % Controls
    insertCD = true;
    
    freq_shift = false;
    
    CDcompensate = true;
    
    if CDcompensate && ~insertCD
        error('Cannot compensate for CD if it is not added in the first place')
    end
    %% 

    % Inserting CD
    D = 17; % Group velocity dispersion in ps/nm/km
    CLambda = 1550e-9; % Central wavelength of optical carrier in m
    linkLength = 40e3; % Fibre link length in m
    
    clear channelOutQAM channelOutQAM_shaped
    
    if insertCD
        channelOutQAM = CDInsertion(modOpticSigQAM, SpS, Rs, D, CLambda, linkLength, 2); 
        channelOutQAM_shaped = CDInsertion(modOpticSigQAM_shaped, SpS, Rs, D, CLambda, linkLength, 2); 
    else
        channelOutQAM = modOpticSigQAM;
        channelOutQAM_shaped= modOpticSigQAM_shaped;
    end
    
    %% 

    % Inserting rotations and/or phase shifts
    theta = 0; phi = 0;
    channelOutQAM = PolRotate(channelOutQAM, theta, phi); % Unshaped
    channelOutQAM_shaped = PolRotate(channelOutQAM_shaped, theta, phi); % Shaped
    %% 

    % Adding AWGN

    SNRb_dB = 8; % SNR per bit in dB
    SNRb_Lin = 10.^(SNRb_dB/10); % linear SNR per symbol

    noisePw = zeros(1,length(lineWidth_Arr)); 
    noisePw_shaped = zeros(1,length(lineWidth_Arr));
    
    channelOutQAM = repelem(channelOutQAM, 1, 1, length(lineWidth_Arr));
    channelOutQAM_shaped = repelem(channelOutQAM_shaped, 1, 1, length(lineWidth_Arr));    
    
    for i=1:length(lineWidth_Arr)
        [channelOutQAM(:,:,i), noisePw(1,i)] = awgn(channelOutQAM(:,:,i), convertSNR(SNRb_dB,'ebno', samplespersymbol=1, bitspersymbol=4), 'measured'); % Unshaped
        [channelOutQAM_shaped(:,:,i), noisePw_shaped(1,i)] = awgn(channelOutQAM_shaped(:,:,i), convertSNR(SNRb_dB,'ebno', samplespersymbol=2, bitspersymbol=4), 'measured'); % Shaped
    end
%% 
    % Frequency Shift
    Elo = Laser(laserPw_dBm, 0, SpS, Rs, NSymb, 2); % LO signal with non-zero linewidth for unshaped signals
    Elo_shaped = Laser(laserPw_dBm, 0, SpS, Rs, samples/SpS, 2); % LO signal with non-zero linewidth for shaped signals
    
    freq_offset = 0;
    if freq_shift && freq_offset == 0
        error('Frequency offset should be non-zero')
    end
    
    Elo_shifted = LOFrequencyShift(Elo, freq_offset, T); % Frequency-shifted LO signal for unshaped information signal
    Elo_shaped_shifted = LOFrequencyShift(Elo_shaped, freq_offset, T); % Frequency-shifted LO signal for shaped information signal
%% 

    % Generating output from optical front-end
    eta_ph = 1; % Responsivity of photodetectors
    
    % Setting size of third dimension
    OpticFEOut_3d = length(lineWidth_Arr);
    
    OpticFEOutQAM = zeros(size(channelOutQAM,1), 4, OpticFEOut_3d);
    OpticFEOutQAM_shaped = zeros(size(channelOutQAM_shaped,1), 4, OpticFEOut_3d);
    
    for i=1:OpticFEOut_3d
        Elo = Laser(laserPw_dBm, lineWidth_Arr(i)*(1-distr), SpS, Rs, NSymb, 2); % LO for unshaped signals
        Elo_shaped = Laser(laserPw_dBm, lineWidth_Arr(i)*(1-distr), SpS, Rs, samples/SpS, 2); % LO for shaped signals
    
        OpticFEOutQAM(:,:,i) = OpticalFrontEnd(channelOutQAM(:,:,i), Elo, eta_ph, 2); % Unshaped
        OpticFEOutQAM_shaped(:,:,i) = OpticalFrontEnd(channelOutQAM_shaped(:,:,i), Elo_shaped, eta_ph, 2); % Unshaped
    end
    
    % ---------------------------------------------------------------------------------------
    
    % TIA amplification and sampling using ADC
    TIA_Gain = 1; % TIA amplification factor
    
    % Amplification and separation of in-phase and quadrature components of the two polarizations
    % 16-QAM
    OpticFEOutQAM_I = TIA_Gain*OpticFEOutQAM(:,1:2:end,:); OpticFEOutQAM_Q = TIA_Gain*OpticFEOutQAM(:,2:2:end,:); % Unshaped 
    OpticFEOutQAM_shaped_I = TIA_Gain*OpticFEOutQAM_shaped(:,1:2:end,:); OpticFEOutQAM_shaped_Q = TIA_Gain*OpticFEOutQAM_shaped(:,2:2:end,:); % Shaped 
%% 

    % Scaling factor between output of IQM and output of optical front-end
    laserEMag = sqrt(laserPw_Lin/2);
    scaleFactor = sqrt(2)*TIA_Gain*eta_ph*laserEMag*laserEMag;
    % scaleFactor=1;
    
    % Cancelling the scaling factor:
    voltageQAM_I = OpticFEOutQAM_I/scaleFactor; voltageQAM_Q = OpticFEOutQAM_Q/scaleFactor; % Unshaped 
    
    % Shaped
    for i=1:OpticFEOut_3d
        voltageQAM_shaped_I(:,:,i) = upfirdn(OpticFEOutQAM_shaped_I(:,:,i)/scaleFactor, filterResponse, 1, 1);
        voltageQAM_shaped_Q(:,:,i) = upfirdn(OpticFEOutQAM_shaped_Q(:,:,i)/scaleFactor, filterResponse, 1, 1);
    end
    voltageQAM_shaped_filtered = voltageQAM_shaped_I(filterSpan*SpS+1:end-filterSpan*SpS+1,:,:)+1i*voltageQAM_shaped_Q(filterSpan*SpS+1:end-filterSpan*SpS+1,:,:);
%% 
    % Frequency Downconversion
    if freq_shift
        for i=1:OpticFEOut_3d
            [voltageQAM_FreqRec(:,:,i), Delta_f(:,i)] = FrequencyRecovery(voltageQAM_I(:,:,i)+1i*voltageQAM_Q(:,:,i), SpS*Rs);
            [voltageQAM_shaped_FreqRec(:,:,i), Delta_f_shaped(:,i)] = FrequencyRecovery(voltageQAM_shaped_filtered(:,:,i), SpS*Rs);
        end
    end
%% 
    % Applying the overlap and save method for CD compensation

    % CD Compensation
    c = 299792458; % Speed of light in m/s
    N_CD = ceil(6.67*(D*1e-6)*(CLambda^2)*linkLength*(Rs^2)*SpS/(2*pi*c)); % CD Equalizer Length#
    N_CD = 50;
    N_FFT = 2^9;
    
    if freq_shift
        CDcompIn = voltageQAM_FreqRec; CDcompIn_shaped = voltageQAM_shaped_FreqRec;
    else
        CDcompIn = voltageQAM_I+1i*voltageQAM_Q; CDcompIn_shaped = voltageQAM_shaped_filtered;
    end
    
    if CDcompensate
        for i=1:OpticFEOut_3d
            voltageQAM_CDC(:,i) = OverlapSaveCDC(CDcompIn(:,1,i), D, linkLength, CLambda, Rs, 1, SpS, N_FFT, N_CD); % Unshaped
            voltageQAM_shaped_CDC(:,i) = OverlapSaveCDC(CDcompIn_shaped(:,1,i), D, linkLength, CLambda, Rs, 1, SpS, N_FFT, N_CD); % Shaped
        end
    end
    %% 

    % Adaptive Equalization

    % Without channel distortion
    NTaps = 5; NRemove = 5000; runs = 1; N1 = runs*length(sourceSymbolsQAM)/2; Mu_T = 1e-2; Mu_DD = Mu_T; Mu_p=1e-2;
    % trainSymbsQAM = codedSymbolsQAM; % trainSymbsQAM(2:2:end) = conj(trainSymbsQAM(2:2:end));
    trainSymbsQAM = sourceSymbolsQAM;
    
    equalizedSymbsQAM = zeros(NSymb*runs, OpticFEOut_3d); equalizedSymbsQAM_shaped = zeros(NSymb*runs, OpticFEOut_3d);
    
    if CDcompensate
        equalizerIn = voltageQAM_CDC;
        equalizerIn_shaped = voltageQAM_shaped_CDC;
    else
        equalizerIn = reshape(voltageQAM_I(:,1,:)+1i*voltageQAM_Q(:,1,:), NSymb*SpS, OpticFEOut_3d, 1)
        equalizerIn_shaped = reshape(voltageQAM_shaped_filtered(:,1,:), NSymb*SpS, OpticFEOut_3d,1);
    end
    
    for i=1:OpticFEOut_3d
        [equalizedSymbsQAM(:,i), ~, ~, ~, ~, ~, ~] = AdapEqualizerAlamouti(equalizerIn(:,i), trainSymbsQAM, SpS, NTaps, Mu_T, Mu_DD, Mu_p, N1, '16QAM', runs); % Unshaped
        [equalizedSymbsQAM_shaped(:,i), ~, ~, ~, ~, ~, ~] = AdapEqualizerAlamouti(equalizerIn_shaped(:,i), trainSymbsQAM, SpS, NTaps, Mu_T, Mu_DD, Mu_p, N1, '16QAM', runs); % Shaped
    end
    equalizedSymbsQAM = equalizedSymbsQAM(length(trainSymbsQAM)*(runs-1)+1:end,:);
    equalizedSymbsQAM_shaped = equalizedSymbsQAM_shaped(length(trainSymbsQAM)*(runs-1)+1:end,:);
    %% 

    % Decision and Decoding
        
    % Final symbols
%     IntsQAM_final = qamdemod(equalizedSymbsQAM(NRemove+1:end-NTaps,:), 16, UnitAveragePower=true);
%     BitsQAM_final = int2bit(IntsQAM_final, 4);
%     [~, BERQAM_final] = biterr(BitsQAM_final, sourceBitsQAM(NRemove*4+1:end-NTaps*4)); % Unshaped
%     
    IntsQAM_shaped_final = qamdemod(equalizedSymbsQAM_shaped(NRemove+1:end-NTaps,:), 16, UnitAveragePower=true);
    BitsQAM_shaped_final = int2bit(IntsQAM_shaped_final, 4);
    [~, BERQAM_shaped_final] = biterr(BitsQAM_shaped_final, sourceBitsQAM(NRemove*4+1:end-NTaps*4)); % Shaped

end