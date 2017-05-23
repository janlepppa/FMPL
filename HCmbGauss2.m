%A heuristic search algorithm for optimizing the MPL. Two steps: 
%Step 1: Performs Markov blanket discovery for each individual node
%Step 2: Combines the MBs into a graph using the AND- and OR-method
%(Described in MPL paper)
%INPUT:     S - d x d unscaled covariance matrix
%           n - sample size
%           prior - flag whether to use prior or not
%           maxParents - maximum number of nodes that can be put in one
%           Markov blanket during search. 
%OUTPUT:    UG_OR - Adjacency matrix for the graph found by the OR-method
%           UG_AND - Adjacency matrix for the graph found by the OR-method
function [UG_OR, UG_AND, MB_times] = HCmbGauss2(S,n, prior,maxParents)
    TIMETEST = 0;
    
    if TIMETEST
        disp('timetest')
    end

    [~, p] = size(S);
  
    MB_times = zeros(p,1);
    UG=zeros(p);
    MBs = cell(1,p);
    

    if TIMETEST
    % perform Markov blanket sequantially
        for j=1:p  
            t1 = tic;
            [MB, ~] = findMBGauss2(j, n, S, prior, maxParents);
            MB_times(j) = toc(t1);
            %UG(MB, j)=1;
            MBs{j} = MB;


        end
    
    else
        % perform Markov blanket discovery in parallel
        parfor j=1:p  

            [MB, ~] = findMBGauss2(j, n, S, prior, maxParents);
            %UG(MB, j)=1;
            MBs{j} = MB;


        end
    end
    
    % combine the found Markov blankets
    for j=1:p  
   
        MB = MBs{j};
        UG(MB, j)=1;
        
    end
    
    % form consistent UGs
    UG=UG+UG';    
    UG_OR=UG>0;    
    UG_AND=(UG==2);
