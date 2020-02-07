function z = shrinkage(x, kappa)
%     z = max( 0, x - kappa ) - max( 0, -x - kappa );
    z = (abs(x) > abs(kappa)) .* (x  - kappa .* x ./ abs(x));
end
