function [u, omega, rho] = gauss_seidel(A, f, tol)
    % Decompose A into L and U (D is included in L)
    D = diag(diag(A));
    L = tril(A);  % Includes D (lower triangular part + diagonal)
    U = triu(A, 1);  % Upper triangular part without the diagonal
    
    % Initial guess for the solution vector
    u = sparse(size(A, 1), 1);
    
    % Iterative process
    max_iter = 10000;  % Maximum iteration limit for safety
    iter = 0;
    while true
        u_old = u;  % Store the old solution
        
        % Update each component of u in-place
        for i = 1:length(u)
            % Calculate the sums for the Gauss-Seidel update
            sum1 = A(i, 1:i-1) * u(1:i-1);  % Using updated values
            sum2 = A(i, i+1:end) * u_old(i+1:end);  % Using old values
            u(i) = (f(i) - sum1 - sum2) / A(i, i);
        end
        
        % Check for convergence
        if norm(u - u_old, inf) < tol  % Use infinity norm for stricter convergence check
            break;
        end
        
        iter = iter + 1;
        if iter > max_iter
            warning('Gauss-Seidel: Maximum number of iterations reached.');
            break;
        end
    end
    
    % The relaxation parameter for Gauss-Seidel is typically 1
    omega = 1;
    
    % Spectral radius calculation (optional, for diagnostics)
    B = inv(D + tril(A, -1)) * triu(A, 1);
    rho = spectral_radius(B);
end
