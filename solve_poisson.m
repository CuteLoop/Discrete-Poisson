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




    A = build_matrix(N);

   


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

