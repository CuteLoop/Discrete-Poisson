% TEST.M - Test the build_matrix function from solve_poisson.m
% This script tests the build_matrix function to ensure it works correctly
% for small grid sizes, displaying properties and verifying matrix dimensions
% and sparsity.

% Define a small grid size for testing
N = 10;

% Test the negative Laplacian case
is_negative = true;
A = build_matrix(N, is_negative);

% Display the properties of the matrix A
fprintf('Matrix A (negative Laplacian) for N=%d:\n', N);
spy(A); % Visual representation of sparsity pattern
title(sprintf('Sparsity Pattern of Matrix A (negative Laplacian) for N=%d', N));
xlabel('Columns');
ylabel('Rows');

% Check if A has the correct dimensions and type
assert(isequal(size(A), [N^2, N^2]), 'Matrix A has incorrect dimensions.');
assert(issparse(A), 'Matrix A is not sparse.');

% Display a message to indicate success
fprintf('Matrix A (negative Laplacian) successfully built with correct dimensions and sparsity for N=%d.\n', N);

% Test the positive Laplacian case
is_negative = false;
A = build_matrix(N, is_negative);

% Display the properties of the matrix A
fprintf('Matrix A (positive Laplacian) for N=%d:\n', N);
spy(A); % Visual representation of sparsity pattern
title(sprintf('Sparsity Pattern of Matrix A (positive Laplacian) for N=%d', N));
xlabel('Columns');
ylabel('Rows');

% Check if A has the correct dimensions and type
assert(isequal(size(A), [N^2, N^2]), 'Matrix A has incorrect dimensions.');
assert(issparse(A), 'Matrix A is not sparse.');

% Display a message to indicate success
fprintf('Matrix A (positive Laplacian) successfully built with correct dimensions and sparsity for N=%d.\n', N);

% Test the matrix for a 3x3 grid, as shown in the provided example
N_test = 3;

% Negative Laplacian case
is_negative = true;
A_test = build_matrix(N_test, is_negative);

% Expected matrix for N=3 (negative Laplacian)
A_expected = [
    -4  1  0  1  0  0  0  0  0;
     1 -4  1  0  1  0  0  0  0;
     0  1 -4  0  0  1  0  0  0;
     1  0  0 -4  1  0  1  0  0;
     0  1  0  1 -4  1  0  1  0;
     0  0  1  0  1 -4  0  0  1;
     0  0  0  1  0  0 -4  1  0;
     0  0  0  0  1  0  1 -4  1;
     0  0  0  0  0  1  0  1 -4
];

% Verify the generated matrix matches the expected matrix
if ~isequal(full(A_test), A_expected)
    fprintf('Matrix A_test (negative Laplacian) does not match the expected matrix for N=%d.\n', N_test);
    fprintf('Generated matrix A_test:\n');
    disp(full(A_test));
    fprintf('Expected matrix A_expected:\n');
    disp(A_expected);
    fprintf('Difference between A_test and A_expected:\n');
    disp(full(A_test) - A_expected);
    error('Matrix A_test (negative Laplacian) does not match the expected matrix for N=3.');
else
    fprintf('Matrix A (negative Laplacian) for N=%d matches the expected matrix.\n', N_test);
end

% Positive Laplacian case
is_negative = false;
A_test = build_matrix(N_test, is_negative);

% Expected matrix for N=3 (positive Laplacian)
A_expected = [
     4  -1   0  -1   0   0   0   0   0;
    -1   4  -1   0  -1   0   0   0   0;
     0  -1   4   0   0  -1   0   0   0;
    -1   0   0   4  -1   0  -1   0   0;
     0  -1   0  -1   4  -1   0  -1   0;
     0   0  -1   0  -1   4   0   0  -1;
     0   0   0  -1   0   0   4  -1   0;
     0   0   0   0  -1   0  -1   4  -1;
     0   0   0   0   0  -1   0  -1   4
];

% Verify the generated matrix matches the expected matrix
if ~isequal(full(A_test), A_expected)
    fprintf('Matrix A_test (positive Laplacian) does not match the expected matrix for N=%d.\n', N_test);
    fprintf('Generated matrix A_test:\n');
    disp(full(A_test));
    fprintf('Expected matrix A_expected:\n');
    disp(A_expected);
    fprintf('Difference between A_test and A_expected:\n');
    disp(full(A_test) - A_expected);
    error('Matrix A_test (positive Laplacian) does not match the expected matrix for N=3.');
else
    fprintf('Matrix A (positive Laplacian) for N=%d matches the expected matrix.\n', N_test);
end

