% Sample data from Gaussian undirected graph.
% INPUT
%   Amatrix     adjancy-matrix of some UG
%   n           sample size
% OUTPUT
%   DATA        n x d sized data matrix
function [DATA] = sampleGaussianData(Amatrix, n)

% number of variables
d = length(Amatrix);

% omega will be the precision matrix of the data to be generated
omega = zeros(d,d);

% lower bound for absolute value of elements
a = 0.1;

% upper bound
b = 0.9;

% probability of off-diagonal element having minus sign
p = 0.5;

% fill with random numbers, leaving off-diagonal zeros that are present in Amatrix in place
for i = 1:d
    for j=1:i
        if Amatrix(i,j) == 1 || i == j
            r = rand;
            
            if r < p && i ~= j
                m = -1;
            else
                m = 1;
            end
            omega(i,j) = m*(a + (b-a)*rand);
            omega(j,i) = omega(i,j);
        end
    end
end

%omega

% parameter to ensure positive definitness
alpha = abs(min(eig(omega))) + 0.1;

% add stuff to diagonal
omega = omega + alpha*eye(d);

% invert to obtain the covariance matrix
covx = inv(omega);

% assume zero mean
meanx = zeros(d,1);

DATA = mvnrnd(meanx, covx, n);
