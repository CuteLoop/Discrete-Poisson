
function [u, omega, rho] = sor(A, f, omega, tol)
    % Determine omega if not provided
    
    N = sqrt(size(A, 1));
        
    
    omega = 1.82;
    
    % Decompose A into D, L, and U
    D = diag(diag(A));
    L = tril(A, -1);
    U = triu(A, 1);
    
    % Compute the iteration matrix B and modified f
    D_omegaL_inv = inv(D + omega * L);
    B = D_omegaL_inv * ((1 - omega) * D - omega * U);
    c = omega * D_omegaL_inv * f;
    
    % Compute the spectral radius
    rho = spectral_radius(B);
    assert(rho < 1, 'SOR: Spectral radius is not < 1');
    
    % Initialize u as a sparse vector
    u = sparse(size(D, 1), 1);
    
    % Iterative process
    while true
        u_new = B * u + c;
        if norm(u_new - u) < tol
            break;
        end
        u = u_new;
    end
end

