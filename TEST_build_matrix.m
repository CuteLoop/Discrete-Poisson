% TEST_build_matrix.m - Test the build_matrix function for sparse and non-sparse matrices

% Test for various grid sizes
test_sizes = [3, 5, 10];
errors_detected = false;

% Function to check if matrix is sparse or non-sparse and run the corresponding tests
function run_tests(N, A)
    h = 1 / (N + 1);
    scale_factor = 1 / h^2;

    if issparse(A)
        fprintf('Matrix A for N = %d is sparse.\n', N);

        % Check if A has the correct dimensions
        assert(isequal(size(A), [N^2, N^2]), 'Matrix A has incorrect dimensions.');
        
        % Display success message
        fprintf('Matrix A (Sparse) for N = %d successfully built with correct dimensions and sparsity.\n', N);

        % Visualize the sparsity pattern for smaller N
        if N <= 5
            spy(A);
            title(sprintf('Sparsity Pattern of Matrix A (Sparse) for N=%d', N));
            xlabel('Columns');
            ylabel('Rows');
            drawnow;
        end

        % Check if the matrix matches the negative Laplacian (scaled or non-scaled)
        if max(abs(diag(A))) == 4 && min(nonzeros(A)) == -1
            fprintf('Matrix A (Sparse) for N = %d matches the non-scaled negative Laplacian.\n', N);
        elseif max(abs(diag(A))) == 4 * scale_factor && min(nonzeros(A)) == -scale_factor
            fprintf('Matrix A (Sparse) for N = %d matches the scaled negative Laplacian.\n', N);
        else
            fprintf('Matrix A (Sparse) for N = %d does not match expected scaled or non-scaled negative Laplacian.\n', N);
            error('Generated sparse matrix does not match the expected type (scaled or non-scaled negative Laplacian).');
        end

    else
        fprintf('Matrix A for N = %d is not sparse (full matrix).\n', N);

        % Check if A has the correct dimensions
        assert(isequal(size(A), [N^2, N^2]), 'Matrix A has incorrect dimensions.');

        % Display success message
        fprintf('Matrix A (Non-Sparse) for N = %d successfully built with correct dimensions.\n', N);

        % Check if the matrix matches the negative Laplacian (scaled or non-scaled)
        if max(abs(diag(A))) == 4 && min(nonzeros(A)) == -1
            fprintf('Matrix A (Non-Sparse) for N = %d matches the non-scaled negative Laplacian.\n', N);
        elseif max(abs(diag(A))) == 4 * scale_factor && min(nonzeros(A)) == -scale_factor
            fprintf('Matrix A (Non-Sparse) for N = %d matches the scaled negative Laplacian.\n', N);
        else
            fprintf('Matrix A (Non-Sparse) for N = %d does not match expected scaled or non-scaled negative Laplacian.\n', N);
            error('Generated non-sparse matrix does not match the expected type (scaled or non-scaled negative Laplacian).');
        end
    end
end

% Run the tests for each grid size
fprintf('\n---- Running Tests for Sparse and Non-Sparse Matrices ----\n');
for N = test_sizes
    try
        % Generate the matrix
        A = build_matrix(N);

        % Run tests based on whether the matrix is sparse or not
        run_tests(N, A);

    catch ME
        fprintf('Error encountered in test for N = %d:\n', N);
        disp(ME.message);
        errors_detected = true;
    end
end

% Final report
if errors_detected
    fprintf('\nTests completed with errors.\n');
else
    fprintf('\nAll tests completed successfully.\n');
end
