function [L U] = factor(C, rho, D)
    L = chol( C'*C + rho*D'*D, 'lower' );
    L = sparse(L);
    U = sparse(L');
end