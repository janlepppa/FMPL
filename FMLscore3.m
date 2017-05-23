% Computes a score based on logarithm of marginal data density defined in formula (24) of Consonni's
% paper. Includes the precomputed sparsity promoting prior term (priorVec)
% INPUT
%  j       index of the current node
%  pa      indices of a proposed parent set of the node j
%  S       p x p - matrix containing sums of squares and products of observations
%           vectors' coordinates i.e. a sample covariance matrix without
%           constant (1/(n-1))
%  priorVec precomputed vector containing prior terms corresponding to number of
%           parents of node. Zero vector implies no prior.
% OUTPUT    
%  score
function [score] =  FMLscore3(j, pa, n, S, priorVec)
% family of node j
fa = [j pa];

% amount of observations and number of variables
[p, ~] = size(S);

% parameter, value p-1 used in Consonni paper
alambda = p-1;

%training sample size, can be set to its minimal value 1, if alambda = p-1
n0 = 1;

% number of parents
pj = length(pa);
 
% Matrix S_papa in the paper, submatrix of S containing terms related to
% variables in vector pa
Spa = S(pa,pa);
 
% % Matrix S_fafa in the paper, submatrix of S containing terms related to
% % variables in vector fa
Sfa = S(fa,fa);

% Compute the score in parts. Pi term
score1 = ( -(n - n0)/2 ) *log(pi);

% (n0/n) fraction part
score2 = ((alambda + n0 - p + 2*pj + 1)/2)*log(n0/n);

% Part involving products of gamma functions. When calculating marginal likelihood, terms cancel out 
% and only the part with j = 1 on the numerator of (24) will remain.
arg1 = (alambda +n -p +pj + 1)/2;
arg2 = (alambda +n0 -p + pj + 1)/2;

score3  = gammaln(arg1) - gammaln(arg2);

% part with determinants of matrix S's submatrices
Lfa = chol(Sfa,'lower');
Lpa = chol(Spa,'lower');

l_detSfa = 2*sum(log(diag(Lfa)));
l_detSpa = 2*sum(log(diag(Lpa)));

score4 = ( -(n-n0)/2 )* (l_detSfa - l_detSpa);

%%%%%%%%%%
% prior term, these are zeros if prior is not used
if pj == 0
    prior = 0;
else
    prior = priorVec(pj);
end

score = prior + score1 + score2 + score3 + score4;
    
end
 

