%%%
% A demo of MATLAB implementation of Algorithm 1 from BCM17
% (as contained in FA.R)
% Written by Martin S. Copenhaver (www.mit.edu/~mcopen)
%%%


%%
% Design parameters

p = 20;
R = 2;

%%

%%%
% Set random set for reproducability, define SDP solver, and turn off
% eigenvalues warnings.
%%%

rng(1,'twister');

cvx_solver Mosek; % other options include SDPT3 and SeDuMi

if true %%% remove various eigenvalue warnings
    warning('off','MATLAB:nargchk:deprecated');
    warning('off','MATLAB:eigs:TooManyRequestedEigsForRealSym');
    warning('off','MATLAB:eigs:TooManyRequestedEigsForComplexNonsym');
end


%% Create true underlying matrices

L = randn(p,R);
TRUE_Theta = L*(L');
Phi = rand(p,1);
TRUE_Phi = trace(TRUE_Theta)/sum(Phi)*Phi;
S = TRUE_Theta+ diag(TRUE_Phi);

%% Test FA functionality

[Theta, Phi] = FA(S,2,1e-5,5000); % first argument is covariance matrix;
                                  % second argument is desired number of
                                  % factors; third argument is numerical
                                  % tolerance for convergence (relative
                                  % optimality gap); and final argument is
                                  % maximum number of iterations.

% distance between Phi and TRUE_Phi

norm(Phi-TRUE_Phi)

% distance between Theta and TRUE_Theta

norm(Theta-TRUE_Theta,'fro')