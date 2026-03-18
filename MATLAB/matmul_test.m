function C = matmult_compare(A, B)
%MATMULT_COMPARE Multiply two square matrices using loops and using A*B.
%   C = MATMULT_COMPARE(A,B) checks inputs, computes product using explicit
%   triple for-loops and using MATLAB built-in multiplication, verifies
%   they match (within floating-point tolerance), and returns the result
%   from the built-in multiplication.
%
%   Inputs:
%     A, B - square numeric matrices of the same size.
%
%   Output:
%     C - product A*B computed with MATLAB's built-in operator.
%
%   The function also prints the time taken by each method and the max
%   absolute difference between results.

% Validate inputs
validateattributes(A, {'numeric'}, {'2d','square'}, mfilename, 'A', 1);
validateattributes(B, {'numeric'}, {'2d','square'}, mfilename, 'B', 2);
if ~isequal(size(A), size(B))
    error('A and B must be the same size.');
end

n = size(A,1);

% Multiply using explicit triple for-loops
C_loop_wrong = zeros(n, n, class(A));
C_loop       = zeros(n, n, class(A));

tic;
for i = 1:n
    for k = 1:n
        for j = 1:n
            C_loop_wrong(i,j) = C_loop_wrong(i,j) + A(i,k) * B(k,j);
        end
    end
end
time_loop_wrong = toc;

tic;
for j = 1:n
    for k = 1:n
        for i = 1:n
            C_loop(i,j) = C_loop(i,j) + A(i,k) * B(k,j);
        end
    end
end
time_loop = toc;

% Multiply using built-in operator
tic;
C = A * B;
time_builtin = toc;

% Compare results
diff_max = max(abs(C(:) - C_loop(:)));
diff_max = max([diff_max, max(abs(C(:) - C_loop_wrong(:)))]);

fprintf('Size: %dx%d\n', n, n);
fprintf('Time (loops wrong order):    %.6f s\n', time_loop_wrong);
fprintf('Time (loops)            :    %.6f s\n', time_loop);
fprintf('Time (builtin)          :    %.6f s\n', time_builtin);
fprintf('Max abs difference      :    %.3e\n'  , diff_max);

% If difference is larger than tolerance, warn
tol = 100 * eps(class(A)) * max(1, max(abs(C(:))));
if diff_max > tol
    warning('Results differ by more than tolerance (%.3e).', tol);
end

end

%% --------------------------------------------------------------------- %%
clear
close all
clc

n = 1024;

% Generate random square matrices A and B
A = rand(n);
B = rand(n);

% Call the matrix multiplication comparison function
C_result = matmult_compare(A, B);