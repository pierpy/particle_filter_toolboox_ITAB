function [ val ] = evalpdfobs( y, y_, noise )
% evalpdfobs Evaluate the pdf from prior observation distribution
    for i = 1:length(y)
        v(i) = pdf(noise.type_y, y(i) - y_(i), 0, noise.sigma_y(i)) + 1e-300;
    end
    val = mean(v);
end

