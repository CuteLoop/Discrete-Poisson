
## Developer Roadmap for Discrete Poisson Solver Project

### 1. **Understanding the Problem**
- **Objective**: Solve the Poisson equation \(\Delta u(x, y) = f(x, y)\) with Dirichlet boundary conditions on the unit square \([0, 1] \times [0, 1]\).
- **Discretization**:
  - Use finite-difference approximation for the Laplacian operator.
  - Implement \( u_{i, j} \) as the solution on a grid \( N \times N \) with spacing \( h = \frac{1}{N+1} \).
  - Ensure boundary conditions \( u(0, y) = u(1, y) = u(x, 0) = u(x, 1) = 0 \) are applied.
- **Matrix Form**:
  - Construct a sparse coefficient matrix \( A \) for the 2D discretized Laplacian using `build_matrix.m`.
  - Use iterative solvers (Jacobi, Gauss-Seidel, SOR) to solve the system \( A u = f \).

### 2. **Function Overview**
- **`solve_poisson.m`**:
  - Main function that interfaces with the iterative solvers.
  - Handles input parsing and calls `build_matrix`, `jacobi`, `gauss_seidel`, and `sor` functions based on user-specified parameters.
- **`build_matrix.m`**:
  - Constructs the sparse matrix \( A \) representing the discretized 2D Laplacian.
- **`TEST.m`**:
  - Runs tests to verify the implementation of the functions and the correctness of the solution.
  - Includes checks for convergence and correct output.

### 3. **Completion Steps**
1. **Implement `build_matrix.m`**:
   - Verify that the matrix \( A \) is constructed using sparse matrices for efficient storage and computation.
   - Ensure that the Kronecker product accurately represents the 2D grid Laplacian structure.
2. **Complete `gauss_seidel.m` and `sor.m`**:
   - Implement the iterative logic for Gauss-Seidel and SOR methods.
   - Ensure that convergence criteria based on the spectral radius are met.
   - Optimize \( \omega \) for the SOR method if not provided.

### 4. **Testing and Validation**
- **Run `TEST.m`**:
  - Verify if the functions return the correct output format and meet convergence criteria.
  - Ensure that the solutions for different grid sizes \( 3 \leq N \leq 40 \) are accurate.
- **Debug with Visualizations**:
  - Use `script.m` to generate plots and validate the solution visually.
  - Check the correctness by comparing the gradient and solution behavior for various methods.

### 5. **Performance Considerations**
- **Sparse Matrices**:
  - Ensure that all matrices are sparse to optimize memory usage and computation.
- **Spectral Radius**:
  - Validate the spectral radius using the `spectral_radius` function to confirm convergence.
- **Optimization**:
  - Implement `optimize_omega` if needed to find the best \( \omega \) for the SOR method.

### 6. **Future Enhancements**
- Add more test cases in `TEST.m` for boundary conditions and larger grid sizes.
- Improve performance by parallelizing parts of the code (e.g., using `parfor` in MATLAB).
- Extend to solve non-uniform grid problems or different boundary conditions.

