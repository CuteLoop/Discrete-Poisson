function A = build_matrix(N, is_negative)
    % Create diagonal and off-diagonal vectors for the tridiagonal blocks
    e = ones(N, 1);
    if is_negative
        T = spdiags([e -4*e e], -1:1, N, N);  % Tridiagonal part for negative Laplacian
        S = spdiags([e e], [-1 1], N, N);  % Off-diagonal part representing neighbors in the other direction
    else
        T = spdiags([-e 4*e -e], -1:1, N, N);  % Tridiagonal part for positive Laplacian
        S = spdiags([-e -e], [-1 1], N, N);  % Off-diagonal part representing neighbors in the other direction
    end
    I = speye(N);  % Identity matrix of size N

    % Construct the full block matrix using Kronecker products
    A = kron(I, T) + kron(S, I);
end
