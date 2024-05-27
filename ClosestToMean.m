function [y] = ClosestToMean(x, NOutput)

    y = zeros(NOutput, size(x,2));

    for i=1:size(x,2)
        indices_w_zero = x(:,i);
        indices_wo_zero = find(indices_w_zero);
        x_wo_zeros = indices_w_zero(indices_wo_zero);
        
        if length(x_wo_zeros) <= NOutput
            y(:,i) = x_wo_zeros;
        else
            x_wo_zeros = sort(x_wo_zeros);
            y(:,i) = x_wo_zeros(1:NOutput,:);
        end
    end
end