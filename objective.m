function p = objective(A, b, lambda, x, D)
    p = ( 1/2*sum(abs(A*x - b).^2) + lambda*norm(D*x,1) );
end