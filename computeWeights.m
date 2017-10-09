function [ ww ] = computeWeights( y, yp, R )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    for i = 1: size(yp, 2)
        ww(i) = mvnpdf(y - yp(:,i), 0,  R) + 1e-300;
    end

end

