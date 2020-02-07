function [s, history] = lasso(C, i, D, lambda, rho, alpha)   
% lasso  Solve lasso problem via ADMM
%
% [s, history] = lasso(C, i, D, lambda, rho, alpha)   ;
%
% Solves the following problem via ADMM:
%
%   minimize 1/2*|| Cr - i ||_2^2 + \lambda ||Dx ||_1
%
% The solution is returned in the vector r.
%
% history is a structure that contains the objective value, the primal and
% dual residual norms, and the tolerances for the primal and dual residual
% norms at each iteration.
%
% rho is the augmented Lagrangian parameter.
%
% alpha is the over-relaxation parameter (typical values for alpha are
% between 1.0 and 1.8).


t_start = tic;

QUIET    = 1;
MAX_ITER = 10000;
ABSTOL   = 1e-4;
RELTOL   = 1e-2;

[m, n] = size(C);

% save a matrix-vector multiply

Atb = C' * i;
r = zeros(n,1);
s = zeros(n,1);
u = zeros(n,1);

% cache the factorization
[L U] = factor(C, rho , D);

if ~QUIET
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
      'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
end

for k = 1:MAX_ITER

    % r-update
    q = Atb + rho*D'*(s - u);  

       r = U \ (L \ q);
    % s-update with relaxation
    sold = s;
    x_hat = alpha*D*r + (1 - alpha)*sold;
    s = shrinkage(x_hat + u, lambda/rho);

    % u-update
    u = u + (x_hat - s);

    % diagnostics, reporting, termination checks
    
    history.objval(k)  = objective(C, i, lambda, r,D);

    history.r_norm(k)  = norm(D*r - s);
    history.s_norm(k)  = norm(-rho*D'*(s - sold));

    history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(D*r), norm(-s));
    history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*D'*u);

    if ~QUIET
        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
            history.r_norm(k), history.eps_pri(k), ...
            history.s_norm(k), history.eps_dual(k), history.objval(k));
    end

    if (history.r_norm(k) < history.eps_pri(k) && ...
       history.s_norm(k) < history.eps_dual(k))
         break;
    end

end

if ~QUIET
    toc(t_start);
end
end

