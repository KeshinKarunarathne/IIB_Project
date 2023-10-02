function [rOut] = GSOP(rIn, NPol)
    
    % Taking the in-phase component as reference:
    rIOrt = rIn(:,1)/sqrt(mean(rIn(:,1).^2));

    % Orthogonalization:
    rQInt = rIn(:,2)-mean(rIn(:,1).*rIn(:,2))*rIn(:,1)/mean(rIn(:,1).^2);
    rQOrt = rQInt/sqrt(mean(rQInt.^2));

    % Complex output signal:
    rOut = [rIOrt rQOrt];

    if NPol == 2
        % Taking the in-phase component as reference:
        rIOrt = rIn(:,3)/sqrt(mean(rIn(:,3).^2));

        % Orthogonalization:
        rQInt = rIn(:,4)-mean(rIn(:,3).*rIn(:,4))*rIn(:,3)/mean(rIn(:,3).^2);
        rQOrt = rQInt/sqrt(mean(rQInt.^2));

        % Complex output signal:
        rOut = [rOut rIOrt rQOrt];
    end
end