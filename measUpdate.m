function [ yp ] = measUpdate(k, predictedParticles, model )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    for i = 1: size(predictedParticles, 2)
        if isempty(model.u)
            yp(:, i) = model.g(k, model.p0, predictedParticles(:, i), []);
        else
            yp(:, i) = model.g(k, model.p0, predictedParticles(:, i), model.u(:,k));
        end
    end
end

