%% Inexact hierarchical scale separation.
%
% Execute 
%
%   p = publish('ihss.m', struct('evalCode', false)); 
%   open(p)
% 
% to show formated documentation.
%
% This MATLAB implementation of IHSS is not performance-optimized and 
% should instead provide an easy way to modify the algorithm for future
% developers.
% 
% The function 
%
%   ihss(A, B, X0, Nel, Nloc, idx_bar, idx_hat, opt)
% 
% approximates the solution $\vec{x}^\ast$ of a linear system $\mathbf{A}\vec{x} = \vec{b}$
% stemming from a discontinuous Galerkin method using modal basis functions such that
%
% $$\frac{\|\mathbf{A}\vec{x}^\ast - \vec{b}\|}{\|\mathbf{A}\vec{x}_0 -\vec{b}\|} < \eta$$
%
% where $\vec{x}_0$ is an initial guess.  For details, see
%
% <html>
%   <p style="border:2px; border-style:solid; border-color:#000000; padding: 1em;">
%     C Thiele, M Araya-Polo, FO Alpak, B Rivière, F Frank<br>
%     Inexact hierarchical scale separation: A two-scale approach for linear
%     systems from discontinuous Galerkin discretizations<br>
%     Computers and Mathematics with Applications, 74(8), 1769–1778, 2017<br>
%     <a href="http://dx.doi.org/10.1016/j.camwa.2017.06.025">DOI: 10.1016/j.camwa.2017.06.025</a>
%     </p>
% </html>
%
% Authors: Florian Frank and Christopher Thiele.
% Contact: florian.frank@fau.de
%
%% Input arguments
% # |A|       - System matrix of size |[Nel*Nloc, Nel*Nloc]|.
% # |B|       - System right-hand side of size |[Nel*Nloc, 1]|.
% # |X0|      - Initial guess of size |[Nel*Nloc, 1]|.
% # |Nel|     - Number of elements.
% # |Nloc|    - Number of local degrees of freedom (assumed to be equal on every element).
% # |idx_bar| - Indices of |X| that correspond to the coarse-scale solution, vector of size |[Nel]|.
% # |idx_hat| - Indices of |X| that correspond to the fine-scale solution, vector of size |[Nel*(Nloc-1)]|.
%
%% Optional input arguments
% Optional arguments are specified by an 8th input argument of type |struct|
% (i.e. they have no order).
% If one of the below fields are not specified, the default value is used.
%
% * |opt.eta|                 - Relative tolerance for the global residual (default: |1E-6|).
% * |opt.nu|                  - Number of fine-scale solves per HSS cycle (default: 8).
% * |opt.m_anderson|          - Sequence length within Anderson acceleration (default: 3).
% * |opt.max_HSS_cycles|      - Maximum allowed number of HSS cycle (default: |Nel|).
% * |opt.CSS_type|            - Name of the coarse-scale solver (|'gmres'| (default), |'cg'|, |'bicgstab'|).
% * |opt.max_CSS_iter|        - Maximum allowed number of coarse-scale solver iterations (default: |Nel|).
% * |opt.coarse_rel_tol_init| - Relative tolerance for the coarse-scale solver in the first HSS cycle (default: 0.1).
% * |opt.is_print|            - Print the global residual norms for every HSS cycle (default: true).
% * |opt.is_plot|             - Plot the residual norms against iteration levels (default: true).
function X = ihss(A, B, X0, Nel, Nloc, idx_bar, idx_hat, ...               % Mandatory arguments.
                  opt)                                                     % Arguments that have a default value.

%% Assertions for mandatory input arguments.
assert(Nel >= 1)
assert(Nloc >= 2)
assert(isequal(size(A),  [Nel*Nloc, Nel*Nloc]))
assert(isequal(size(B),  [Nel*Nloc, 1]))
assert(isequal(size(X0), [Nel*Nloc, 1]))
assert(length(idx_bar) == Nel)
assert(length(idx_hat) == Nel*(Nloc - 1))
assert(nargin == 7 || nargin == 8, 'ihss has seven mandatory and one optional input argument, see ''doc ihss''.')

%% Setting unspecified optional input parameters (i.e. completing struct 'opt').
if nargin == 8
  assert(isstruct(opt), 'The argument for optional parameters must be a struct.')
  if ~isfield(opt, 'eta')
    opt.eta = 1E-6;
  end
  if ~isfield(opt, 'nu')
    opt.nu = 8;
  end
  if ~isfield(opt, 'm_anderson')
    opt.m_anderson = 3;
  end
  if ~isfield(opt, 'max_HSS_cycles')
    opt.max_HSS_cycles = Nel;
  end
  if ~isfield(opt, 'CSS_type')
    opt.CSS_type = 'gmres';
  end
  if ~isfield(opt, 'max_CSS_iter')
    opt.max_CSS_iter = Nel;
  end
  if ~isfield(opt, 'coarse_rel_tol_init')
    opt.coarse_rel_tol_init = 0.1;
  end
  if ~isfield(opt, 'is_print')
    opt.is_print = true;
  end
  if ~isfield(opt, 'is_plot')
    opt.is_plot = true;
  end
end

%% Assertions for optional arguments.
assert(opt.coarse_rel_tol_init > 0)
assert(opt.eta > 0)
assert(opt.nu >= 1)
assert(opt.m_anderson >= 0)

%% Extracting the system blocks.
% In performance-relevant implementations, this should be avoided by having
% the data ready in the desired format.
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
                       norm(C_hat*X_bar_iter + A_hat*X_hat_old - B_hat)^2);% Global residual of initial guess.
global_res_cur = global_res_init;                                          % Current global residual.
if opt.is_print
  fprintf('Residual %f\n', global_res_cur);
end
num_CSS_iter_list    = [];
global_res_norm_list = [];                                                 % Arrays to store the residual norms to plot them later.
coarse_res_norm_list = [];
fine_res_norm_list   = [];

while (global_res_cur/global_res_init >= opt.eta)                          % Check global relative residual norm.
  if k > opt.max_HSS_cycles                                                % Check if we exceed the number of allowed HSS cycles.
    error('IHSS did not converge within the allowed number of cycles.')
  end
  
  %% Step 2.1: Initialization.
  k = k + 1;
  if opt.is_print
    fprintf('HSS cycle number %d.\n', k);
  end

  %% Step 2.6: Determine a tolerance for the coarse-scale solver.
  % Given A*x = b, MATLAB gmres provides only a stopping criterion for a
  % relative residual norm of the form ||A*x_k - b||/||b|| < tol.  Hence,
  % we need to convert the absolute tolerance for the coarse-scale solver 
  % to the above one.
  if k == 1 
    coarse_rel_tol = opt.coarse_rel_tol_init;
  else
    coarse_rel_tol = fine_res_norm/norm(B_bar - C_bar*X_hat_old);
  end
  
  %% Step 2.2: Coarse system approximation.
  % Solve the coarse scale system with the coarse-scale solver (CSS).
  % If solver converges print information, otherwise terminate.
  switch opt.CSS_type
    case 'gmres'
      gmres_restart = 30;
      [X_bar_iter, flag, relres, num_CSS_iter] = ...
        gmres(A_bar, B_bar - C_bar*X_hat_old, gmres_restart, coarse_rel_tol, opt.max_CSS_iter, [], [], X_bar_iter);
      num_CSS_iter = num_CSS_iter(2);
    case 'cg'
      [X_bar_iter, flag, relres, num_CSS_iter] = ...
        pcg(A_bar, B_bar - C_bar*X_hat_old, coarse_rel_tol, opt.max_CSS_iter, [], [], X_bar_iter);
    case 'bicgstab'
      [X_bar_iter, flag, relres, num_CSS_iter] = ...
        bicgstab(A_bar, B_bar - C_bar*X_hat_old, coarse_rel_tol, opt.max_CSS_iter, [], [], X_bar_iter);  
    otherwise
      error('The value for opt.CSS_type is unknown.')
  end % switch
  % Store number of CSS iterations for plot.
  num_CSS_iter_list = [num_CSS_iter_list, num_CSS_iter];                   %#ok<AGROW>
  if (flag == 0)
    if opt.is_print
      fprintf('  CSS (%s) iterations: %.1f (rel. res.: %.3e).\n', opt.CSS_type, num_CSS_iter, relres);
    end
  else
    error('CSS (%s) diverged (error flag: %d)', opt.CSS_type, flag)
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
  [~, X_hat_new]  = evalc('AndAcc(g_fun, X_hat_old, opt.m_anderson, opt.nu, atol, rtol);'); 
  coarse_res_norm = norm(A_bar*X_bar_iter + C_bar*X_hat_new - B_bar);
  fine_res_norm   = norm(C_hat*X_bar_iter + A_hat*X_hat_new - B_hat);
  global_res_cur  = sqrt(coarse_res_norm^2 + fine_res_norm^2);
  % Store residual norms for plotting.
  global_res_norm_list = [global_res_norm_list, global_res_cur];           %#ok<AGROW>
  coarse_res_norm_list = [coarse_res_norm_list, coarse_res_norm];          %#ok<AGROW>
  fine_res_norm_list   = [fine_res_norm_list,   fine_res_norm];            %#ok<AGROW>
  if opt.is_print
    fprintf('  Abs. residual norm: %.3e\n', global_res_cur);
    fprintf('  Rel. residual norm: %.3e\n', global_res_cur/global_res_init);
  end
  X_hat_old = X_hat_new;
end % while

% Plot the residual norms for the fine-scale, coarse-scale, and original system.
if opt.is_plot
  % Left y-axis for logarithmic residuals.
  yyaxis left
  semilogy(global_res_norm_list, 'k');
  hold on,  box on,  grid on
  semilogy(coarse_res_norm_list, 'g--');
  semilogy(fine_res_norm_list, 'r--');
  % Right axis is number of CSS iterations per HSS cycle.
  yyaxis right
  plot(num_CSS_iter_list)
  ylim([0 ceil(max(num_CSS_iter_list))])
  legend('global residual norm', 'coarse-scale residual norm', 'fine-scale residual norm', 'CSS iterations', ...
         'Location', 'southwest')
  xlabel('HSS cycle number')
end 

%% Reordering of the solution vector X.
% In performance-relevant implementations, this should be avoided by having
% the data ready in the desired format.
X = zeros(Nel*Nloc, 1);
X(idx_bar) = X_bar_iter;
X(idx_hat) = X_hat_new;

end % function
