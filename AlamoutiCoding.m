function [codedSymbols] = AlamoutiCoding(inputSymbols)
    % This function implements Alamouti encoding for a single stream of
    % symbols

    codedSymbols = zeros(length(inputSymbols), 2);

    for i=1:2:length(inputSymbols)
        codedSymbols(i,:) = [inputSymbols(i,1) inputSymbols(i+1,1)];
        codedSymbols(i+1,:) = [-1*conj(inputSymbols(i+1,1)) conj(inputSymbols(i,1))];
    end

end
