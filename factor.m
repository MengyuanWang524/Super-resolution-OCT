function [L U] = factor(A, rho, D)
    L = chol( A'*A + rho*D'*D, 'lower' );
    % force matlab to recognize the upper / lower triangular structure
    L = sparse(L);
    U = sparse(L');
end