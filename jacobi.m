
function [u, omega, rho] = jacobi(A, f, tol)
    D = diag(diag(A));
    L = tril(A,-1);
    U = triu(A,1);    
    Dinv = inv(D);
    B = Dinv*(U+L);
    f = Dinv*f;
    rho = spectral_radius(B);
    assert(rho < 1, 'Jacobi: Spectral radius is not < 1');
    u = sparse(size(D,1),1);
    while true
        u_new = f - B * u;
        if norm(u_new - u) < tol
            break;
        end
        u = u_new;
    end
    omega = NaN;
end