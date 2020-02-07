function p = objective(C, i, lambda, r, D)
    p = ( 1/2*sum(abs(C*r - i).^2) + lambda*norm(D*r,1) );
end