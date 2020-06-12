function s = shrinkage(r, kappa)
%     s = max( 0, r - kappa ) - max( 0, -r - kappa );
    s = (abs(r) > abs(kappa)) .* (r  - kappa .* r ./ abs(r));
end
