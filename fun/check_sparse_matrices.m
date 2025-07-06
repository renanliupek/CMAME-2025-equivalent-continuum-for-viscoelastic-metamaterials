function check_sparse_matrices(M, C, K)
    % Function to check if M is Hermitian positive definite
    % and if C and K are Hermitian positive semidefinite.
    % This version is optimized for large and sparse matrices.
    
    fprintf('Checking matrix M:\n');
    check_sparse_hermitian_positive_definite(M);
    
    fprintf('\nChecking matrix C:\n');
    check_sparse_hermitian_positive_semidefinite(C);
    
    fprintf('\nChecking matrix K:\n');
    check_sparse_hermitian_positive_semidefinite(K);
end

function check_sparse_hermitian_positive_definite(A)
    % Check if A is Hermitian
    if issymmetric(A)
        fprintf('Matrix is Hermitian.\n');
    else
        fprintf('Matrix is NOT Hermitian.\n');
        return; % If not Hermitian, no need to check further
    end
    
    % Check if A is positive definite using Cholesky decomposition
    [R, p] = chol(A);
    if p == 0
        fprintf('Matrix is positive definite.\n');
    else
        fprintf('Matrix is NOT positive definite. Cholesky decomposition failed.\n');
    end
end

function check_sparse_hermitian_positive_semidefinite(A)
    % Check if A is Hermitian
    if issymmetric(A)
        fprintf('Matrix is Hermitian.\n');
    else
        fprintf('Matrix is NOT Hermitian.\n');
        return; % If not Hermitian, no need to check further
    end
    
    % Try using eigs with a small diagonal shift to improve convergence
    opts.tol = 1e-8; % Relax tolerance
    opts.maxit = 1000; % Increase iterations
    shift = 1e-5; % Small shift to aid convergence
    
    try
        % Try shifted matrix
        min_eig = eigs(A + shift * speye(size(A)), 1, 'smallestreal', opts);
        if min_eig - shift >= 0
            fprintf('Matrix is positive semidefinite (with small shift).\n');
        else
            fprintf('Matrix is NOT positive semidefinite. Smallest eigenvalue (shifted): %f\n', min_eig - shift);
        end
    catch ME
        fprintf('Failed to compute eigenvalues using shifted eigs. Reason: %s\n', ME.message);
        % Fallback to full matrix eig if small enough
        if nnz(A) < 1e6
            try
                full_min_eig = min(eig(full(A)));
                if full_min_eig >= 0
                    fprintf('Matrix is positive semidefinite (full matrix eigenvalues).\n');
                else
                    fprintf('Matrix is NOT positive semidefinite. Smallest eigenvalue: %f\n', full_min_eig);
                end
            catch ME2
                fprintf('Full matrix eigenvalue check also failed. Reason: %s\n', ME2.message);
            end
        else
            fprintf('Matrix is too large for full eigenvalue computation.\n');
        end
    end
end
