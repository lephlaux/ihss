format compact                                                             % Omit empty lines in output.
close all                                                                  % Close all figures.

Nel = 32^3;                                                                % Number of elements.
Nloc = 4;                                                                  % Number of local degrees of freedom.
idx_bar = 1 : Nloc : Nel*Nloc;                                             % Indices for the averages.                                             
idx_hat = setdiff(1 : Nel*Nloc, idx_bar);                                  % Indices for higher oder terms.
eta = 1E-10;                                                               % Relative tolerance.
nu = 8;                                                                    % Fine-scale solves per HSS cycle.

load('system_good.mat', 'A', 'B')                                          % Droplet scenario, small time step.                           
% load('system_bad.mat', 'A', 'B')                                         % Spinodal decomposition, large time step.
X0 = zeros(size(B));                                                       % Initial guess.

X_HSS = ihss(A, B, X0, Nel, Nloc, idx_bar, idx_hat, eta, nu);              % Call solver.
