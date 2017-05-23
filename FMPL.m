% Computes OR, AND and HC (if on) graphs.
% S - unscaled covariane matrix of data
% n - sample size
% prior - 1/0, whether to use prior
% HCON - 1/0,  whether to compute HC
% maxParents - maximum number of nodes that can be put in one
%           Markov blanket during search. 
% OUTPUT: undirected graphs
function [OR, AND, HC, OAtime, HCtime, MBtimes] = FMPL(S, n, prior, HCON, maxParents)

[~, d] = size(S);

% if prior is not on, (log) prior terms will be zeros
priorVec = zeros(d,1);

%vector containing possible numbers of parents/sizes of Markov
%blankets
pjs = 1:d-1;

% precompute prior values
if prior == 1    
    ms = (pjs + 1).*(pjs)/2;
    
    %a,b parameters for a beta distribution (vectors that repeat these values)
    a = 0.5 * ones(size(pjs));
    b = 0.5 * ones(size(pjs));
    
    % vector containing prior penalties
    priorVec = betaln(a+pjs,b+ms-pjs) - betaln(a,b);
    
end
% compute OR and AND graphs
t2 = tic;
[OR, AND, MBtimes] = HCmbGauss2(S, n, priorVec,maxParents);
OAtime = toc(t2);
% ...  and maybe HC
if HCON == 1
    
    t3 = tic;
    HC = ORHcGauss2(OR,S, n, priorVec);
    HCtime = toc(t3);
else
    HC = [];
end
