function [ predictedParticles ] = predictParticles(k, pastParticles, model)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    w = model.Q*randn(size(pastParticles,1), size(pastParticles,2));
    for i = 1: size(pastParticles, 2)
        if isempty(model.u)
            predictedParticles(:, i) = pastParticles(:, i) + (model.f(k, model.p0, pastParticles(:, i), []) + w(:, i)) * model.Dt;
        else
            predictedParticles(:, i) = pastParticles(:, i) + (model.f(k, model.p0, pastParticles(:, i), model.u(:,k)) + w(:, i)) * model.Dt;
        end
    end

end

