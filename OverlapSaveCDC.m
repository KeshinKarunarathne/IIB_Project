function [Out] = OverlapSaveCDC(In,D,L,CLambda,Rs,NPol,SpSIn,NFFT,NOverlap)
    
    % Parameters:
    c = 299792458 ;D=D*1e-6;

    % Index for coefficient calculation and Nyquist frequency:
    n = (-NFFT/2:NFFT/2-1)' ; fN = SpSIn*Rs/2;

    % Calculating the CD frequency response:
    HCD = exp(-1i*pi*CLambda^2*D*L/c*(n*2*fN/NFFT).^2);

    if NPol == 2
        HCD = cat(3,HCD,HCD);
    end

    % Guaranteeing that NOverlap is even:
    NOverlap = NOverlap + mod(NOverlap,2);

    % Extending the input signal so that the blocks can be properly formed:
    AuxLen = size(In,1)/(NFFT-NOverlap);
    
    if AuxLen ~= ceil(AuxLen)
        NExtra = ceil(AuxLen)*(NFFT-NOverlap)-size(In,1);
        In = [In(end-NExtra/2+1:end,:); In ; In(1:NExtra/2,:)];
    else
        NExtra = NOverlap;
        In = [In(end-NExtra/2+1:end,:); In ; In(1:NExtra/2,:)];
    end

    % Blocks:
    BlocksV = reshape(In(:,1),NFFT-NOverlap,size(In,1)/(NFFT-NOverlap));

    if NPol == 2
        BlocksH = reshape(In(:,2),NFFT-NOverlap,size(In,1)/(NFFT-NOverlap));
        Blocks = cat(3,BlocksV,BlocksH) ; clearvars BlocksV BlocksH In;
    else
        Blocks = BlocksV ; clearvars BlocksV In;
    end

    % Preallocating the output blocks and the overlap for the first block:
    Out = zeros(size(Blocks)) ; Overlap = zeros(NOverlap,1,NPol);

    % Compensating for the chromatic dispersion:
    for i = 1:size(Blocks,2)
        % Input block with overlap:
        InB = [Overlap ; Blocks(:,i,:)];

        % FFT of the input block:
        InBFreq = fftshift(fft(InB));

        % Filtering in frequency domain:
        OutFDEFreq = InBFreq.*HCD;

        % IFFT of the block after filtering:
        OutFDE = ifft(ifftshift(OutFDEFreq));

        % Overlap:
        Overlap = InB(end-NOverlap+1:end,1,:);

        % Output block:
        OutB = OutFDE(NOverlap/2+1:end-NOverlap/2,1,:);

        % Assigning the samples to the output signal:
        Out(:,i,:) = OutB;
    end

    % Output column vector:
    OutV = reshape(Out(:,:,1),numel(Out(:,:,1)),1);
    if NPol == 2
        OutH = reshape(Out(:,:,2),numel(Out(:,:,2)),1);
        Out = [OutV OutH];
    else
        Out = OutV;
    end

    % Quantity of samples to discard:
    DInit = 1+(NExtra+NOverlap)/2 ; DFin = (NExtra-NOverlap)/2;
    
    % Removing the overlapped zeros:
    Out = Out(DInit:end-DFin,:);
end





