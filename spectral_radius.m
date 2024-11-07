function rho = spectral_radius(A)
    rho = abs( eigs(A,1,'largestabs','MaxIterations',512) );
end