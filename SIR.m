function U = SIR(phi,x0,sub)

%%%%%%%%%%%%%%%%
%% Parameters %%
%%%%%%%%%%%%%%%%
if sub==0
    % Default values
    Rfac=0.5;      % Reduction of R at each iteration
    dPdx=0.95;     % Initial values of R
else
    % Values when subiterating
    Rfac=0.8;      % Reduction of R at each iteration
    dPdx=0.9999;   % Initial values of R
end

N=numel(x0);     % Number of equations
tol=1e-8;        % Solution accuracy
imax=100;         % Maximum number of iterations
Js=1000;           % Maximum number of subiterations
a_c=2;           % Critical magnitude for alpha
S1_min=-5e-2;    % Parameter for monotonicity check

%%%%%%%%%%%%%%%%%%%%
%% Initialization %%
%%%%%%%%%%%%%%%%%%%%

x0 = x0(:);             % Turn x0 into column vector
xold = zeros(size(x0)); % Store the old x0
R = dPdx*ones(N,1);     % R0

J = @(x)jacobian(phi,x,N);

%%%%%%%%%%%%%%%%%%%%
%% Iteration Loop %%
%%%%%%%%%%%%%%%%%%%%
for n = 1:imax

    % --------------------
    % I: Compute new alpha
    % --------------------

    x    = x0;            % Initial guess
    phi0 = (phi(x))';     % Evaluate initial guess
    J0   = J(x);          % Evaluate center diff if no Jacobian
    
    % Test for singularity
    % Convert to full matrix, in case J is sparse.
    condJ = rcond(full(J0));  
    if condJ <= eps
        error('Jacobian singular in SIR');
    end

    % Compute R-I
    RI0 = RI_proc(R,N);
    % Compute alpha
    alpha = alpha_proc(RI0,J0,N);
    % --------------
    % II: Compute x1
    % --------------
    x1 = Phi_proc(x0,phi0,alpha);

    %%%%%%%%%%%%%%%%%%%
    %% Subiterations %%
    %%%%%%%%%%%%%%%%%%%
    % Only subiterate if step length is increasing
    if sub && max(abs(double(x1-x0)) - abs(double(x0-xold))) > 0
        mon = x0-x1;

        % -----------------------
        % V: Check validity of x1
        % -----------------------
        for j=1:Js
            % Evaluate phi(x1)
            x = x1;            % initial guess
            phi1 = (phi(x))';  % evaluate initial guess

            S1 = mon .* (x1-Phi_proc(x1,phi1,alpha));
            S2 = max(abs(alpha),[],2);
           
            % Perform all comparisons at once, 
            % store in boolean vector ifR
            % A particular row of ifR will be true if subiteration is
            % required in that dimension.
            ifR = (S2 >= a_c | S1 < S1_min);

            % Break if no subiterations are needed
            if max(ifR) == 0
                break;
            end

            % Increase R towards I
            for i=1:N
                % Increase i:th dimension
                if ifR(i) == 1
                    R(i) = (3*R(i)+1)*0.25;
                end
            end
            % Update RI
            RI0 = RI_proc(R,N);

            % Recompute alpha
            alpha = alpha_proc(RI0,J0,N);
            x1 = Phi_proc(x0,phi0,alpha);
        end
        % -----------------------
        % Subiterations end here:
        % -----------------------
    end

    % ------------------------
    % Accuracy test and update
    % ------------------------
    eps_acc = mean(abs(x1-x0));
    xold = x0; % Backup x0 so that it can be used to compare step lengths.
    x0 = x1;
    if eps_acc < tol
        break
    end

    % ----------
    % Decrease R
    % ----------
    % This was moved to the end of the iteration step, rather than being
    % performed at the beginning of all steps but the first. This removes
    % the need to test if n>1 at each step.
    R = Rfac*R;
end
U = x0';

%%%%%%%%%%%%%%%
%% Functions %%
%%%%%%%%%%%%%%%
function alpha = alpha_proc(RI,J,N)
alpha = eye(N)+RI/J;

function Phi = Phi_proc(x,phi,alpha)
% Computes Phi (x and phi are allowed to be matrices)
Phi = alpha*(x-phi)+phi;

function RI = RI_proc(R,N)
% Computes RI
RI = zeros(N,N);
for i=1:N
    RI(i,i) = R(i)-1;
end

function J = jacobian(phi,x,N)
h=eps^(1/3);       % Finite difference length
Nf=numel(phi(x));  % Degrees of freedrom
J=zeros(N,Nf);     % Allocate J
for i=1:N
    dx=zeros(N,1); dx(i)=dx(i)+h;
    J(:,i)=(phi(x+dx)-phi(x-dx))/h/2; % Center difference
end

% If J is a cell array, we extract its component parts
if iscell(J)
    r = double(J{1});
	c = double(J{2});
	J = J{3};
	sparse_def = true; % Flag sparse version
else
	sparse_def = false; % Flag dense version
end
 
% Transform J0 into I-J0
if sparse_def
	% Sparse version
	J = sparse(r,c,J,N,N);
	J = speye(N)-J;
else
    % Dense version
	J = eye(N)-J;
end


