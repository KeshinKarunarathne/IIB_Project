function [y,pArr,w11_Arr,w12_Arr,w21_Arr,w22_Arr,errorArr,BER] = ExpDataProcessing(x, filter, Rs, trainSymbs, D, CLambda, linkLength, N_CD, N_FFT, NTaps, Mu_T, Mu_DD, Mu_p, N1, modFormat, NRemoveF, NRemoveB, runs)

    % This function processes data collected from experiments

    % Inputs
    % x - Input symbols as a single-column vector at 256 GSa/s
    % filter - RRC filter for matched filtering
    % Rs - Symbol rate in Hz
    % trainSymbs - original source symbols as a single coulmn vector
    % D - Fibre dispersion constant in ps/nm/km
    % CLambda - Central wavelength of laser
    % linkLength - length of fibre in m
    % N_CD - Length of overlap for overla-and-save
    % N_FFT - FFT size for overlap-and-save
    % NTaps - Number of taps per filter
    % Mu_T - step-size for the training stage
    % Mu_DD - step-size for the decision-directed stage
    % Mu_p - step-size for phase correction
    % N1 - Number of iterations of the training stage before switching to decision-directed mode
    % modFormat - Modulation format
    % NRemoveF - Symbols not to consider at the front due to errors when calculating BER
    % NRemoveB - Symbols not to consider at the end due to errors when calculating BER

    % Outputs
    % y - Equalized symbols as a single-column vector
    % pArr - variation of p as a single-column vector
    % w11_Arr - variation of coefficient of w11 
    % w12_Arr - variation of coefficient of w12 
    % w21_Arr - variation of coefficient of w21 
    % w22_Arr - variation of coefficient of w22 
    % errorArr - variation of error as a two-column matrix
    % BER - bit error rate of equalized symbols

    % ---------------------------------------------------------------------

    % Donwconversion to baseband using circshift
    delta = 10; % to cut out the DC component from analysis
    spec = fft(x); % spectrum of signal at IF
    [~,index] = max(abs(spec(delta+1:end-delta)));
    shift  = delta + index;
    shifted_spec = circshift(spec,shift); % spectrum centred around 0 Hz
    baseband_sig = ifft(shifted_spec); % baseband signal

    % ---------------------------------------------------------------------

    % Resampling from 256GSa/s to 100GSa/s
    resampled_sig = resample(baseband_sig, 100, 256);

    % ---------------------------------------------------------------------

    % Matched Filtering
    filtered_sig = upfirdn(resampled_sig, filter, 1, 1);

    % ---------------------------------------------------------------------

    % Frame Synchronisation
    SpS = 2;

    % Finding samples to extract using cross-correlation
    [c,lags] = xcorr(filtered_sig, repelem(trainSymbs, SpS, 1));
    
    synch_index = lags(c == max(c)); % Lag with largest cross-correlation
    
    % Extracting synchronised samples
    synchronised_sig = filtered_sig(synch_index:synch_index + SpS*length(trainSymbs)-1,:);
    
    % Normalising the extracted samples
    norm_sig = synchronised_sig/max(abs(synchronised_sig));

    % ---------------------------------------------------------------------

    % CD compensation
    signal_CDC = OverlapSaveCDC(norm_sig, D, linkLength, CLambda, Rs, 1, SpS, N_FFT, N_CD);

    % ---------------------------------------------------------------------

    % Adaptive Equalization
    % Best results with reruns obtained for norm_sig and normalised symbs
    [y, pArr, w11_Arr,w12_Arr,w21_Arr,w22_Arr,errorArr] = AdapEqualizerAlamouti(signal_CDC, trainSymbs, SpS, NTaps, Mu_T, Mu_DD, Mu_p, N1, modFormat, runs); % Using CDC
%     [y, pArr, w11_Arr,w12_Arr,w21_Arr,w22_Arr,errorArr] = AdapEqualizerAlamouti(norm_sig, trainSymbs/sqrt(10), SpS, NTaps, Mu_T, Mu_DD, Mu_p, N1, modFormat, runs); % Without CDC
    y = y(length(trainSymbs)*(runs-1)+1:end,:);

    % ---------------------------------------------------------------------

    % Decision and Decoding (for 16-QAM)
     
    sourceInts = qamdemod(trainSymbs(NRemoveF+1:end-NRemoveB,1), 16);
    sourceBits = int2bit(sourceInts, 4);
    
    EqualizedInts = qamdemod(y(NRemoveF+1:end-NRemoveB,1), 16); % Change this based on whether s_qam or s_qam/sqrt(10)
    EqualizedBits = int2bit(EqualizedInts, 4);
    
    [~, BER] = biterr(sourceBits, EqualizedBits);

end