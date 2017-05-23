%Discovers the Markov blanket of given node.
%INPUT:     NODE - node index
%           n - sample size
%           S - covariance matrix of data without scaling
%           prior - 1 use prior, 0 no
%           maxParents - max number of elements allowed to be included in
%                       Markov blanket during search
%OUTPUT:    MB_TOP - indices of Markov blanket discovered by the algorithm
%           MB_TOP_SCORE - score of discovered Markov blanket
function [MB_TOP,MB_TOP_SCORE] = findMBGauss2(NODE,n,S, prior, maxParents)



%d = size(DATA,2);
[~, d] = size(S);

% candidates for the Markov blanket (node under consideration removed)
MB_POT=int16(1:d)';
MB_POT(NODE,:)=[];

% initialize the score using an empty Markov blanket
MB_TOP_SCORE = FMLscore3(NODE, [], n, S, prior);
MB_TOP=[];
CONT=1;

while CONT                                                              %adding a node to the MB if it results in an improvement of the score
    CONT=0;
    
    % vector to collect the scores
    MB_CAND_SCORES=zeros(size(MB_POT,1),1);
    
    %Create Markov blanket candidates by adding one node
    %to the existing one.
    %Condition added to control the maximum MB size!!!
    if size(MB_TOP,2)>0 && length(MB_TOP)+1 <= maxParents
        MB_CAND=[MB_TOP(ones(size(MB_POT,1),1),:) MB_POT];
        
        % Max number of elements in MB reached, terminate, no need to check if
        %removing any increases the score since this was done in last iteration
        %elseif length(MB_TOP)+1 >= maxParents  !!!!
    elseif length(MB_TOP)+1 > maxParents
        break;
        
    else
        MB_CAND=MB_POT;
    end
    
    %         if size(MB_TOP,2)>0 && length(MB_TOP)+1 <= maxParents
    %             MB_CAND=[MB_TOP(ones(size(MB_POT,1),1),:) MB_POT];
    %         else
    %             MB_CAND=MB_POT;
    %         end
    
    % score the candidates
    for j=1:size(MB_POT,1)
        MB_CAND_SCORES(j,1) = FMLscore3(NODE, MB_CAND(j,:), n, S, prior);
    end
    
    % score the candidate blankets
    [MB_CAND_TOPSCORE,MB_CAND_TOPLOC]=max(MB_CAND_SCORES);
    
    % if some of the new blankets increases the score, take it as
    % as the new Markov blanket
    if MB_CAND_TOPSCORE > MB_TOP_SCORE
        MB_TOP=MB_CAND(MB_CAND_TOPLOC,:);
        MB_POT(MB_CAND_TOPLOC,:)=[];
        MB_TOP_SCORE=MB_CAND_TOPSCORE;
        CONT=1;
    end
    
    MB_SIZE=size(MB_TOP,2);
    if MB_SIZE>2 && CONT==1                                             %removing a node from the MB if it results in an improvement of the score
        DEL=1;
        while DEL
            DEL=0;
            MB_C=MB_TOP;
            MB_CAND_SCORES=zeros(MB_SIZE,1);
            for j=1:MB_SIZE
                MB_CAND=MB_C;
                MB_CAND(:,j)=[];
                MB_CAND_SCORES(j,1) = FMLscore3(NODE, MB_CAND, n, S, prior);
            end
            [MB_CAND_TOPSCORE,MB_CAND_TOPLOC]=max(MB_CAND_SCORES);
            if MB_CAND_TOPSCORE > MB_TOP_SCORE
                MB_TOP(:,MB_CAND_TOPLOC)=[];
                MB_SIZE=size(MB_TOP,2);
                MB_TOP_SCORE=MB_CAND_TOPSCORE;
                if MB_SIZE>2
                    DEL=1;
                end
            end
        end
    end
end

end

