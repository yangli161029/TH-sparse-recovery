function y = shrink(x, t)
%SHRINK  Soft-thresholding (proximal operator of the L1 norm).
%
%   y = SHRINK(x, t) returns sign(x) .* max(|x| - t, 0), the proximal
%   operator of t * ||.||_1 evaluated at x.
%
%   x and t may be vectors/matrices of compatible sizes.

    y = sign(x) .* max(abs(x) - t, 0);
end
