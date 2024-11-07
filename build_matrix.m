 function A = build_matrix(N)
    % Define constants based on grid spacing
    h = 1 / (N + 1);
    a = 4 / h^2;
    b = -1 / h^2;
    
    % Total number of unknowns in the grid
    M = N^2;
    
    % Initialize sparse matrix
    A = sparse(M, M);
    
    % Loop over each grid point (i, j)
    for i = 1:N
        for j = 1:N
            % Convert (i, j) to linear index
            idx = (i - 1) * N + j;
            
            % Set the main diagonal entry
            A(idx, idx) = a;
            
            % Set the horizontal neighbors if they exist
            if i > 1  % left neighbor
                A(idx, idx - N) = b;
            end
            if i < N  % right neighbor
                A(idx, idx + N) = b;
            end
            
            % Set the vertical neighbors if they exist
            if j > 1  % bottom neighbor
                A(idx, idx - 1) = b;
            end
            if j < N  % top neighbor
                A(idx, idx + 1) = b;
            end
        end
    end
end
