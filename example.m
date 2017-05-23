% A really simple example showing the usage of FMPL function

% the true inverse of covariance
IC = [1 0 0.5; ...
      0 1 0;   ...
      0.5 0 1];

% the true covariance matrix
C = inv(IC);

%sample a data set
n = 100;
D = mvnrnd(zeros(3,1),C,n);

% center the data to have zero mean and scale it to have a standard deviation of one
Dcs = zscore(D);

% unscaled sample covariance matrix
S = Dcs'*Dcs;

% use prior
prior = 1;

% compute HC graph
HCON = 1;

% max number of nodes in a Markov blanket during search
maxParents = n -1; % does not have any effect here since n >> p = 3

%learn the graphs
[OR, AND, HC] = FMPL(S, n, prior, HCON, maxParents)


% true adjacency matrix
trueGraph = (IC -diag(diag(IC))) > 0



