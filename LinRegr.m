function [y_regr] = LinRegr(x, y)
%     slope = x\y;

    ones_vec = [ones(length(x),1) x];
    b = ones_vec\y;

    y_regr = ones_vec*b;
end