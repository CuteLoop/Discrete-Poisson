function best_method = estimate_spectral_radius_comparison(A, max_iter, tol)
    % estimate_spectral_radius_comparison - Estimate and compare the spectral radius using different methods
    % 
    % Syntax: best_method = estimate_spectral_radius_comparison(A, max_iter, tol)
    %
    % Inputs:
    %   A - The matrix for which to estimate the spectral radius
    %   max_iter - Maximum number of iterations (for iterative methods)
    %   tol - Tolerance for convergence (for iterative methods)
    %
    % Outputs:
    %   best_method - The name of the method that gives the closest estimate to the exact value
    
    if nargin < 2
        max_iter = 1000;  % Default maximum iterations
    end
    if nargin < 3
        tol = 1e-6;  % Default tolerance
    end
    
    methods = {'eigs', 'power', 'gershgorin', 'rayleigh', 'norm'};
    rho_values = zeros(length(methods), 1);
    
    % Calculate the spectral radius using each method
    for i = 1:length(methods)
        method = methods{i};
        try
            rho_values(i) = estimate_spectral_radius(A, method, max_iter, tol);
            fprintf('Spectral radius using %s method: %.6f\n', method, rho_values(i));
        catch ME
            fprintf('Error with method %s: %s\n', method, ME.message);
            rho_values(i) = NaN;  % Assign NaN if the method fails
        end
    end
    
    % Find the best estimation (using 'eigs' as the baseline if available)
    rho_eigs = rho_values(strcmp(methods, 'eigs'));
    if ~isnan(rho_eigs)
        % Calculate the absolute difference between 'eigs' and other methods
        differences = abs(rho_values - rho_eigs);
        [~, best_index] = min(differences);
        best_method = methods{best_index};
        fprintf('The method that gives the closest estimation to the eigs method is: %s\n', best_method);
    else
        fprintf('No baseline (eigs) available; cannot determine the best method.\n');
        best_method = 'unknown';
    end
end

function rho = estimate_spectral_radius(A, method, max_iter, tol)
    if nargin < 3
        max_iter = 1000;  % Default maximum iterations
    end
    if nargin < 4
        tol = 1e-6;  % Default tolerance
    end
    
    switch lower(method)
        case 'eigs'
            % Method 1: Use MATLAB's built-in eigs function
            opts = struct('tol', tol, 'maxit', max_iter);
            rho = abs(eigs(A, 1, 'largestabs', opts));
            
        case 'power'
            % Method 2: Power Iteration
            x = rand(size(A, 1), 1);  % Random initial vector
            for k = 1:max_iter
                x_new = A * x;
                rho = norm(x_new);
                x_new = x_new / rho;  % Normalize
                if norm(A * x_new - rho * x_new) < tol
                    break;
                end
                x = x_new;
            end
            
        case 'gershgorin'
            % Method 3: Gershgorin Circle Theorem
            radii = sum(abs(A), 2) - abs(diag(A));
            rho = max(abs(diag(A)) + radii);
            
        case 'rayleigh'
            % Method 4: Rayleigh Quotient Iteration
            x = rand(size(A, 1), 1);  % Random initial vector
            x = x / norm(x);
            for k = 1:max_iter
                mu = (x' * A * x) / (x' * x);  % Rayleigh quotient
                y = (A - mu * eye(size(A))) \ x;  % Solve for the next iteration
                x = y / norm(y);
                if abs((x' * A * x) / (x' * x) - mu) < tol
                    break;
                end
            end
            rho = abs((x' * A * x) / (x' * x));
            
        case 'norm'
            % Method 5: Norm-Based Estimation
            rho = norm(A, 2);  % 2-norm estimate as an upper bound
            
        otherwise
            error('Invalid method specified. Choose from ''eigs'', ''power'', ''gershgorin'', ''rayleigh'', or ''norm''.');
    end
end

A = rand(5);  % Example matrix
best_method = estimate_spectral_radius_comparison(A);
