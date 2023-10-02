function [r] = ADC(i, SPSIn, NPol, ParamFilter, ParamADC)
    
    if ~strcmp(ParamFilter.Type,'NoFilter')
        switch ParamFilter.Type
            case 'RRC'
                % Obtaining the filter transfer function:
                g = RRC(ParamFilter.Span,SpSIn,ParamFilter.Rolloff);

                % Filtering the signals:
                i(:,1)=conv(i(:,1),g,'same'); i(:,2)=conv(i(:,2),g,'same');

                if NPol == 2
                    i(:,3)=conv(i(:,3),g,'same'); i(:,4)=conv(i(:,4),g,'same');
                end
            otherwise
                error('Filter type not supported');
        end
    end

    % Normalizing the signals to unitary energy:
    i(:,1) = i(:,1)/sqrt(mean(abs(i(:,1)).^2));
    i(:,2) = i(:,2)/sqrt(mean(abs(i(:,2)).^2));

    if NPol == 2
        i(:,3) = i(:,3)/sqrt(mean(abs(i(:,3)).^2));
        i(:,4) = i(:,4)/sqrt(mean(abs(i(:,4)).^2));
    end

    % Downsampling the signals to 2 samples per symbol:
    PhaseError = 0; % Default value;
    if isfield(ParamADC,'PhaseError')
        PhaseError = ParamADC.PhaseError*SpSIn;
    end

    Len = length(i(:,1)) ; 
    Pos = PhaseError+(1:(SpSIn/ParamADC.SpS):Len)';

    r(:,1) = interp1(1:Len,i(:,1),Pos,'spline','extrap');
    r(:,2) = interp1(1:Len,i(:,2),Pos,'spline','extrap');

    if NPol == 2
        % Downsampling the signals to 2 samples per symbol:
        r(:,3) = interp1(1:Len,i(:,3),Pos,'spline','extrap');
        r(:,4) = interp1(1:Len,i(:,4),Pos,'spline','extrap');
    end

    if isfield(ParamADC,'FreqError')
        Len = length(r(:,1)) ; Pos = (1:1-ParamADC.FreqError:Len)';
        rAux(:,1) = interp1(1:Len,r(:,1),Pos,'spline','extrap');
        rAux(:,2) = interp1(1:Len,r(:,2),Pos,'spline','extrap');

        if NPol == 2
            rAux(:,3) = interp1(1:Len,r(:,3),Pos,'spline','extrap');
            rAux(:,4) = interp1(1:Len,r(:,4),Pos,'spline','extrap');
        end
        r = rAux;
    end
end
