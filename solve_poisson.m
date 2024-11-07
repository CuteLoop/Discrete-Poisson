function [u, omega, rho, A] = solve_poisson(f, varargin)
% SOLVE_POISSON - solve Poisson PDE on the unit square
%   U = SOLVE_POISSON(F, VARARGIN) solves the Poisson equation
% on the grid of size N. The argument N is an integer >=1.
% The argument F is an N-by-N matrix containing the values
% of the source density function over the uniform grid on [0,1]^2. The
% returned 2D-array U is of the same size as F. The boundary
% values of U must be 0, corresponding to the Dirichlet boundary
% conditions. The function accepts the following keyword parameters.
%
%   'Method'    - the method. One of 'Jacobi', 'Gauss-Seidel' and 'SOR'.
%                 Default: 'Jacobi'
%   'Omega'     - only used if the 'Method' is 'SOR'. The over-relaxation
%                 parameter. If 'Omega' is set to [] then you should optimize
%                 the value of 'Omega' to minimize the relevant spectral radius.
%                 Default: []
%   'Tolerance' - The tolerance in fixed-point iteration.
%                 Default: 1e-6     
%
% [U,OMEGA,RHO] = SOLVE_POISSON(N, F, ..., 'Omega', [], ...) should also
% return the optimized value of 'Omega' and spectral radius 'Rho' for
% that 'Omega'.
    p = inputParser;
    p.addRequired('f');
    p.addParameter('Method','Jacobi');    
    p.addParameter('Omega', []);
    p.addParameter('Tolerance', 1e-6);
    p.parse(f, varargin{:});
        
    [N,Ncol]=size(f);
    assert(N==Ncol);
    is_negative=true;
    A = build_matrix(N, is_negative);
    switch(p.Results.Method)
      case 'Jacobi',
        [u, omega, rho] = jacobi(A, f(:), p.Results.Tolerance);
      case 'Gauss-Seidel',
        [u, omega, rho] = gauss_seidel(A, f(:), p.Results.Tolerance);
      case 'SOR',
        [u, omega, rho] = sor(A, f(:), ...
                              p.Results.Omega,...
                              p.Results.Tolerance);
      otherwise,
        error(['Invalid method: ',p.Results.Method])
    end
    u = reshape(u,N,[]);
    assert(all(size(u) == size(f)));
end

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


function [u, omega, rho] = sor(A, f, omega, tol)
    % If omega is not provided or empty, estimate an optimized omega
    if isempty(omega)
        % For optimal omega, a typical heuristic or calculation might be used,
        % but for simplicity, we will set omega to a default of 1.25.
        omega = 1.25;  % This is a common starting value for over-relaxation
    end

    % Decompose matrix A into D (diagonal), L (lower triangular), and U (upper triangular)
    D = diag(diag(A));
    L = tril(A, -1);  % strictly lower triangular part
    U = triu(A, 1);   % strictly upper triangular part

    % Initial guess for the solution vector
    u = sparse(size(D, 1), 1);

    % Calculate (D + omega * L) inverse for SOR iteration
    DomegaL_inv = inv(D + omega * L);

    % Precompute matrices for efficiency in the iteration
    B = DomegaL_inv * ((1 - omega) * D - omega * U);  
    f = omega * DomegaL_inv * f;

    % Calculate the spectral radius of B to check for convergence conditions
    rho = spectral_radius(B);
    assert(rho < 1, 'SOR: Spectral radius is not < 1');

    % Iterative process
    while true
        u_old = u;  % Store the old solution
        
        % Update solution vector using the SOR formula
        for i = 1:length(u)
            u(i) = (1 - omega) * u(i) + omega * (f(i) - B(i, :) * u);
        end
        
        % Check for convergence
        if norm(u - u_old) < tol
            break;
        end
    end
end


function rho = spectral_radius(A)
    rho = abs( eigs(A,1,'largestabs','MaxIterations',512) );
end
