

function [u, omega, rho] = gauss_seidel(A, f, tol)
    tol = tol/10
    % Decompose matrix A into D (diagonal), L (lower triangular), and U (upper triangular)
    D = diag(diag(A));
    L = tril(A, -1);  % strictly lower triangular part
    U = triu(A, 1);   % strictly upper triangular part

    % Initial guess for the solution vector
    u = sparse(size(D, 1), 1);
    
    % Calculate (D + L) inverse for Gauss-Seidel iteration
    DL_inv = inv(D + L);
    
    % Precompute the matrix B and modified f for efficiency in the iteration
    B = -DL_inv * U;  
    f = DL_inv * f;
    
    % The spectral radius is generally used to check convergence conditions.
    % Calculating the spectral radius of B (optional for Gauss-Seidel)
    rho = spectral_radius(B);
    assert(rho < 1, 'Gauss-Seidel: Spectral radius is not < 1');
    
    % Iterative process
    while true
        u_old = u;  % Store the old solution
        
        % Update solution vector using the Gauss-Seidel formula
        for i = 1:length(u)
            u(i) = f(i) - B(i, :) * u;
        end
        
        % Check for convergence
        if norm(u - u_old) < tol
            break;
        end
    end
    
    % The relaxation parameter for Gauss-Seidel is typically 1
    omega = 1;
end
