function UG_HC = ORHcGauss2(UG_OR,S, n, prior)
%Performs a global greedy hill climb search on the reduced model space created from UG_OR 
%INPUT:     UG_OR - Adjacency matrix for a graph defining the model space
%           S - d x d unscaled covariance matrix
%           n - sample size
%           prior - 1/0, whether to use sparsity promoting prior
%OUTPUT:    UG_HC - Adjacency matrix for the graph discovered by the HC algoritm
    [~, d] = size(UG_OR);
    
    UG_HC=zeros(d);                                                         %the empty graph is set as the initial state
    EDGES_N=sum(UG_OR(:))/2;
    EDGES=zeros(EDGES_N,2);
    pos=1;
    
    for i=1: d                                                              %lists the potential edges according to UG_OR
        MB_TEMP=find(UG_OR(:,i));
        MB_TEMP(MB_TEMP<i,:)=[];
        EDGES(pos:(pos-1+size(MB_TEMP,1)),:)=[i*ones(size(MB_TEMP,1),1) MB_TEMP];
        pos=pos+size(MB_TEMP,1);
    end
    
    EDGE_CHANGE_IMP=zeros(EDGES_N,1);
    EDGE_CHANGE=ones(EDGES_N,1);
    CONT=1;
    
    while CONT                                                             %keeps applying changes (add/remove edge) until no further improvement can be achieved by a local change
        CONT=0;
        EDGE_CHANGE_IND=find(EDGE_CHANGE);
        for i=1:size(EDGE_CHANGE_IND,1)                                     %re-evaluates the nodes affected by the last change    
            EDGE_TEMP=EDGES(EDGE_CHANGE_IND(i),:);
            if UG_HC(EDGE_TEMP(1),EDGE_TEMP(2))==0                          % if edge is not in the graph, add it                           
                NODE=EDGE_TEMP(1);
                MB=find(UG_HC(NODE,:));
                SWOE1=FMLscore3(NODE,MB,n, S, prior);
                MB=[MB EDGE_TEMP(2)];
                SWE1=FMLscore3(NODE,MB,n, S, prior);
                NODE=EDGE_TEMP(2);
                MB=find(UG_HC(NODE,:));
                SWOE2=FMLscore3(NODE,MB,n, S, prior);
                MB=[MB EDGE_TEMP(1)];
                SWE2=FMLscore3(NODE,MB,n, S, prior);
                EDGE_CHANGE_IMP(EDGE_CHANGE_IND(i),1)=SWE1+SWE2-SWOE1-SWOE2;
                EDGE_CHANGE(EDGE_CHANGE_IND(i),1)=0;
            else
                NODE=EDGE_TEMP(1);                                          
                MB=find(UG_HC(NODE,:));
                SWE1=FMLscore3(NODE,MB,n, S, prior);
                MB(:,MB==EDGE_TEMP(2))=[];
                SWOE1=FMLscore3(NODE,MB,n, S, prior);               
                NODE=EDGE_TEMP(2);
                MB=find(UG_HC(NODE,:));
                SWE2=FMLscore3(NODE,MB,n, S, prior);
                MB(:,MB==EDGE_TEMP(1))=[];
                SWOE2=FMLscore3(NODE,MB,n, S, prior);                
                EDGE_CHANGE_IMP(EDGE_CHANGE_IND(i),1)=SWOE1+SWOE2-SWE1-SWE2;
                EDGE_CHANGE(EDGE_CHANGE_IND(i),1)=0;
            end    
        end
        [MAXIMP,MAXIMP_LOC]=max(EDGE_CHANGE_IMP);                           %finds and applies the local change that induces the maximum improvement to the score
        if MAXIMP>0            
            EDGE_TEMP=EDGES(MAXIMP_LOC,:);
            if UG_HC(EDGE_TEMP(1),EDGE_TEMP(2))==0
                UG_HC(EDGE_TEMP(1),EDGE_TEMP(2))=1;
                UG_HC(EDGE_TEMP(2),EDGE_TEMP(1))=1;
            else
                UG_HC(EDGE_TEMP(1),EDGE_TEMP(2))=0;
                UG_HC(EDGE_TEMP(2),EDGE_TEMP(1))=0;
            end
            EDGE_CHANGE(EDGES(:,1)==EDGE_TEMP(1)|EDGES(:,2)==EDGE_TEMP(1),1)=1;
            EDGE_CHANGE(EDGES(:,1)==EDGE_TEMP(2)|EDGES(:,2)==EDGE_TEMP(2),1)=1;
            CONT=1;
        end
    end
%     toc
end
