function [ lambda ] = projectionlambda( lambda )
    if lambda > 1
        lambda = 1;
    elseif lambda < 0
        lambda = 0;
    end
end

