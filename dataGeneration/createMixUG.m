% Creates adjacency matrices of the graphs used in the structure learning tests. Parameter d is the dimension 
% which has to a multiple of 64.  
function [expUG] = createMixUG(d)

loopy = createUG('loopy');
hub = createUG('hub');
grid1 = createUG('grid');
clique = createUG('clique');

UG = blkdiag(loopy, hub, grid1, clique);

expUG = [];

n = d /64;

for i = 1:n
    expUG = blkdiag(expUG, UG);
end
