%% Inexact hierarchical scale separation.
%
% Type 'publish('ihss.m')' to show formated documentation.
%
% This MATLAB implementation of IHSS is not performance-optimized and 
% should instead provide an easy way to modify the algorithm for future
% developers.
% 
% The function 
%
%   'ihss(A, B, X0, Nel, Nloc, idx_bar, idx_hat)'
% 
% approximates the solution $\vec{x}^\ast$ of a linear system $\mathbf{A}\vec{x} = \vec{b}$
% stemming from a discontinuous Galerkin method using modal basis functions such that
%
% $$\frac{\|\mathbf{A}\vec{x}^\ast - \vec{b}\|}{\|\mathbf{A}\vec{x}_0 -\vec{b}\|} < \eta$$
%
% where $\vec{x}_0$ is an initial guess.  For details, see
% _Thiele, Araya-Polo, Alpak, Riviere, Frank - 2017 - Inexact hierarchical scale separation: 
% A two-scale approach for linear systems from discontinuous Galerkin discretizations_.
%
% Contact: florian.frank@rice.edu
%
%% Input arguments
% * |A|       - System matrix of size |[Nel*Nloc, Nel*Nloc]|.
% * |B|       - System right-hand side of size |[Nel*Nloc, 1]|.
% * |X0|      - Initial guess of size |[Nel*Nloc, 1]|.
% * |Nel|     - Number of elements.
% * |Nloc|    - Number of local degrees of freedom (assumed to be equal on every element).
% * |idx_bar| - Indices of |X| that correspond to the coarse-scale solution, vector of size |[Nel]|.
% * |idx_hat| - Indices of |X| that correspond to the fine-scale solution, vector of size |[Nel*(Nloc-1)]|.
%
%% Optional input arguments
% * |eta|                 - Relative tolerance for the global residual (default: |1E-6|).
% * |nu|                  - Number of fine-scale solves per HSS cycle (default: 8).
% * |m_anderson|          - Sequence length within Anderson acceleration (default: 3).
% * |max_HSS_iter|        - Maximum allowed number of HSS cycle (default: |Nel|).
% * |coarse_rel_tol_init| - Relative tolerance for the coarse-scale solver in the first HSS cycle (default: 0.1).
function X = ihss(A, B, X0, Nel, Nloc, idx_bar, idx_hat, ...               % Mandatory arguments.
                  eta, nu, m_anderson, max_HSS_iter, coarse_rel_tol_init)  % Arguments that have a default value.

%% Assertions for mandatory input arguments.
assert(Nel >= 1)
assert(Nloc >= 2)
assert(isequal(size(A),  [Nel*Nloc, Nel*Nloc]))
assert(isequal(size(B),  [Nel*Nloc, 1]))
assert(isequal(size(X0), [Nel*Nloc, 1]))
assert(length(idx_bar) == Nel)
assert(length(idx_hat) == Nel*(Nloc - 1))

%% Setting unspecified optional input parameters.
if nargin < 7,  error('ihss has seven mandatory input arguments, see ''doc ihss''.'); end
if nargin < 8,  eta = 1E-6; end
if nargin < 9,  nu  = 8; end
if nargin < 10, m_anderson = 3; end
if nargin < 11, max_HSS_iter = Nel; end
if nargin < 12, coarse_rel_tol_init = 0.1; end

%% Assertions for optional arguments.
assert(coarse_rel_tol_init > 0)
assert(eta > 0)
assert(nu >= 1)
assert(m_anderson >= 0)

%% Extracting the system blocks (should be done in-situ in a C++ implementation).
A_bar = A(idx_bar, idx_bar);                                               % Extraction of submatrices.
C_bar = A(idx_bar, idx_hat);
C_hat = A(idx_hat, idx_bar);
A_hat = A(idx_hat, idx_hat);
kk = mat2cell(1 : size(A_hat, 1), 1, ones(1, Nel)*(Nloc - 1));             % Get blockdiagonal of A_hat.
tt = cellfun(@(x) A_hat(x, x), kk, 'un', 0);                               % (taken from www.mathworks.com/matlabcentral/answers/39467-extract-diagonal-blocks-of-sparse-matrix).
A_hat_diag = blkdiag(tt{:});
A_hat_off = A_hat - A_hat_diag;
B_bar = B(idx_bar);
B_hat = B(idx_hat);

%% Initialization of the IHSS scheme.
X_bar_iter = X0(idx_bar);
X_hat_old  = X0(idx_hat);
k = 0;
global_res_init = sqrt(norm(A_bar*X_bar_iter + C_bar*X_hat_old - B_bar)^2 + ...
                       norm(C_hat*X_bar_iter + A_hat*X_hat_old - B_hat)^2);
global_res_cur = global_res_init;
fprintf('Residual %f\n', global_res_cur);
global_res_norm_list = [];                                                 % Arrays to store the residual norms to plot them later.
coarse_res_norm_list = [];
fine_res_norm_list   = [];

while (global_res_cur/global_res_init >= eta)                              % Check global relative residual norm.
  if k > max_HSS_iter                                                      % Check if we exceed the number of allowed HSS cycles.
    error('IHSS did not converge within the allowed number of cycles.')
  end
  
  %% Step 2.1: Initialization.
  k = k + 1;
  fprintf('Iteration number: %d.\n', k);

  %% Step 2.6: Determine a tolerance for the coarse-scale solver.
  % Given A*x = b, MATLAB gmres provides only a stopping criterion for a
  % relative residual norm of the form ||A*x_k - b||/||b|| < tol.  Hence,
  % we need to convert the absolute tolerance for the coarse-scale solver 
  % to the above one.
  if k == 1 
    coarse_rel_tol = coarse_rel_tol_init;
  else
    coarse_rel_tol = fine_res_norm/norm(B_bar - C_bar*X_hat_old);
  end
  
  %% Step 2.2: Coarse system approximation.
  % Solve the coarse scale system with (restarted) GMRES.  If solver 
  % converges print information, otherwise terminate.
  gmres_restart = 30;
  [X_bar_iter, flag, relres, numiter] = ...
    gmres(A_bar, B_bar - C_bar*X_hat_old, gmres_restart, coarse_rel_tol, 1000, [], [], X_bar_iter);
  if (flag == 0)
    fprintf('  GMRES(%d) iterations: %d (rel. res.: %e).\n', gmres_restart, numiter(2), relres);
  else
    error('Coarse-scale solver diverged.')
  end
  
  %% Steps 2.3, 2.4: Fine system approximation by Anderson acceleration.
  % 
  % For this step an implementation of the Anderson acceleration is
  % required.  We refer to the implementation of Homer F. Walker in
  %
  %   https://users.wpi.edu/~walker/Papers/anderson_accn_algs_imps.pdf
  %
  % where a Matlab function 'AndAcc' is provided.  We apply the backslash
  % opertator to the global fine-scale system that appears in the fixed-
  % point map 'g_fun'.  In a performance-relevant implementation of IHSS,
  % we recommend solving each of the (small), block-systems by
  % QR-decomposition.  The Anderson acceleration reduces to the (standard)
  % fixed-point iteration if 'm_anderson = 0' making nu iterations.
  d     = B_hat - C_hat*X_bar_iter;
  g_fun = @(X_hat_old) A_hat_diag\(d - A_hat_off*X_hat_old);               %#ok<NASGU>
  atol  = -inf; rtol = -inf;                                               %#ok<NASGU>
  % We use 'evalc' to suppress the output of the Anderson acceleration function.
  [~, X_hat_new]  = evalc('AndAcc(g_fun, X_hat_old, m_anderson, nu, atol, rtol);'); 
  coarse_res_norm = norm(A_bar*X_bar_iter + C_bar*X_hat_new - B_bar);
  fine_res_norm   = norm(C_hat*X_bar_iter + A_hat*X_hat_new - B_hat);
  global_res_cur  = sqrt(coarse_res_norm^2 + fine_res_norm^2);
  % Store residual norms for plotting.
  global_res_norm_list = [global_res_norm_list, global_res_cur];           %#ok<AGROW>
  coarse_res_norm_list = [coarse_res_norm_list, coarse_res_norm];          %#ok<AGROW>
  fine_res_norm_list   = [fine_res_norm_list,   fine_res_norm];            %#ok<AGROW>
  fprintf('  Abs. residual norm: %.3e\n', global_res_cur);
  fprintf('  Rel. residual norm: %.3e\n', global_res_cur/global_res_init);
  X_hat_old = X_hat_new;
end % while

% Plot the residual norms for the fine-scale, coarse-scale, and original system.
semilogy(global_res_norm_list, 'b');
hold on,  box on,  grid on
semilogy(coarse_res_norm_list, 'g--');
semilogy(fine_res_norm_list, 'r--');
legend('global residual norm', 'coarse-scale residual norm', 'fine-scale residual norm');

%% Reordering of the solution vector X.
X = zeros(Nel*Nloc, 1);
X(idx_bar) = X_bar_iter;
X(idx_hat) = X_hat_new;

end % function
