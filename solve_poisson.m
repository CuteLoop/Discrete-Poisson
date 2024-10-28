% SOLVE_POISSON - solve Poisson PDE on the unit square
%   U = SOLVE_POISSON(F, VARARGIN) solves the Poisson equation
% on the grid of size N. The argument F is an N-by-N matrix containing the values
% of the source density function over the uniform grid on [0,1]^2. The
% returned 2D-array U is of the same size as F. The boundary
% values of U must be 0, corresponding to the Dirichlet boundary
% conditions.

function [u, omega, rho, A] = solve_poisson(f, varargin)
    p = inputParser;
    p.addRequired('f');
    p.addParameter('Method', 'Jacobi');
    p.addParameter('Omega', []);
    p.addParameter('Tolerance', 1e-6);
    p.addParameter('MaxIterations', 5000);
    p.parse(f, varargin{:});

    [N, Ncol] = size(f);
    assert(N == Ncol, 'Input f must be a square matrix.');
    
    % Compute grid spacing h
    h = 1 / (N + 1);
    
    % Build the coefficient matrix A with scaling 
    is_negative = false;
    A = build_matrix(N, is_negative) * (1 / h^2);

    switch (p.Results.Method)
        case 'Jacobi'
            [u, omega, rho] = jacobi(A, f(:), p.Results.Tolerance, p.Results.MaxIterations);
        case 'Gauss-Seidel'
            [u, omega, rho] = gauss_seidel(A, f(:), p.Results.Tolerance, p.Results.MaxIterations);
        case 'SOR'
            [u, omega, rho] = sor(A, f(:), p.Results.Omega, p.Results.Tolerance, p.Results.MaxIterations);
        otherwise
            error(['Invalid method: ', p.Results.Method])
    end

    u = reshape(u, N, N);
    assert(all(size(u) == size(f)));
end

% Jacobi method
function [u, omega, rho] = jacobi(A, f, tol, max_iters)
    D = diag(diag(A));
    L = tril(A, -1);
    U = triu(A, 1);
    Dinv = inv(D);
    B = -Dinv * (L + U);
    f_new = Dinv * f;
    rho = spectral_radius(B);
    assert(rho < 1, 'Jacobi: Spectral radius is not < 1');
    u = zeros(size(D, 1), 1);
    iter_count = 0;

    while true
        u_new = B * u + f_new;
        if norm(u_new - u, inf) < tol || iter_count >= max_iters
            break;
        end
        u = u_new;
        iter_count = iter_count + 1;
    end
    omega = NaN;
end

% Gauss-Seidel method
function [u, omega, rho] = gauss_seidel(A, f, tol, max_iters)
    L = tril(A, -1);
    U = triu(A, 1);
    D = diag(diag(A));
    LDinv = inv(L + D);
    B = -LDinv * U;
    f_new = LDinv * f;
    rho = spectral_radius(B);
    assert(rho < 1, 'Gauss-Seidel: Spectral radius is not < 1');
    u = zeros(size(D, 1), 1);
    iter_count = 0;

    while true
        u_new = B * u + f_new;
        if norm(u_new - u, inf) < tol || iter_count >= max_iters
            break;
        end
        u = u_new;
        iter_count = iter_count + 1;
    end
    omega = 1;
end

% SOR method
function [u, omega, rho] = sor(A, f, omega, tol, max_iters)
    if isempty(omega)
        % Perform a linear search for the optimal omega
        best_omega = 1.0;
        best_iter_count = Inf;
        omega_range = 1.0:0.05:1.95;
        for test_omega = omega_range
            [~, iter_count, ~] = sor_iteration(A, f, test_omega, tol, 1000);
            if iter_count < best_iter_count
                best_iter_count = iter_count;
                best_omega = test_omega;
            end
        end
        omega = best_omega;
        fprintf('Optimal omega found: %.2f\n', omega);
    end

    % Now solve using the optimal omega
    [u, ~, rho] = sor_iteration(A, f, omega, tol, max_iters);
end

% Helper function for SOR iterations
function [u, iter_count, rho] = sor_iteration(A, f, omega, tol, max_iters)
    D = diag(diag(A));
    L = tril(A, -1);
    U = triu(A, 1);
    M = (D / omega) + L;
    N_mat = ((1 - omega) / omega) * D - U;
    M_inv = inv(M);
    B = -M_inv * N_mat;
    f_new = M_inv * f;
    rho = spectral_radius(B);

    % Initialize u with zeros for consistency
    u = zeros(size(D, 1), 1);
    iter_count = 0;

    while true
        u_new = B * u + f_new;
        current_error = norm(u_new - u, inf);

        if current_error < tol || iter_count >= max_iters
            break;
        end
        u = u_new;
        iter_count = iter_count + 1;
    end
end

% Spectral radius calculation
function rho = spectral_radius(B)
    rho = max(abs(eigs(B, 1, 'largestabs', 'MaxIterations', 1000)));
end
