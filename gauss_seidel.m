function [u, omega, rho] = gauss_seidel(A, f, tol)
    % Decompose A into D, L, and U
    D = diag(diag(A));
    L = tril(A, -1);
    U = triu(A, 1);
    
    % Compute the iteration matrix B and modified f
    DL_inv = inv(D + L);
    B = DL_inv * U;
    c = DL_inv * f;
    
    % Compute the spectral radius
    rho = spectral_radius(B);
    assert(rho < 1, 'Gauss-Seidel: Spectral radius is not < 1');
    
    % Initialize u as a sparse vector
    u = sparse(size(D, 1), 1);
    
    % Iterative process
    while true
        u_new = c - B * u;
        if norm(u_new - u) < tol
            break;
        end
        u = u_new;
    end
    
    omega = NaN; % Not used in Gauss-Seidel
end