function UG = createUG(type)
% Creates UGs (d=16) in simulation section of MPL article.
% INPUT: type - type of graphs: 'grid','hub','loopy','clique'
% OUTPUT: UG
    if strcmp(type,'grid')
        % 4x4 grid network
        UG=zeros(16);
        UG(1,[2 5])=1;
        UG(2,[1 6 3])=1;
        UG(3,[2 7 4])=1;
        UG(4,[3 8])=1;
        UG(5,[1 6 9])=1;
        UG(6,[2 5 7 10])=1;
        UG(7,[3 6 8 11])=1;
        UG(8,[4 7 12])=1;
        UG(9,[5 10 13])=1;
        UG(10,[6 9 11 14])=1;
        UG(11,[7 10 12 15])=1;
        UG(12,[8 11 16])=1;
        UG(13,[9 14])=1;
        UG(14,[10 13 15])=1;
        UG(15,[11 14 16])=1;
        UG(16,[12 15])=1;
        UG=(UG+UG')>0;
    elseif strcmp(type,'hub')
        % double hub node network
        UG=zeros(16);
        UG(1,2:9)=1;
        UG(9,10:16)=1;
        UG=(UG+UG')>0;
    elseif strcmp(type,'loopy')
        % loopy network
        UG=zeros(16);
        UG(1,[2 3 5 10])=1;
        UG(2,[1 3 6 16])=1;
        UG(3,[1 2 4 9])=1;
        UG(4,[3 5])=1;
        UG(5,[4 1])=1;
        UG(6,[2 7])=1;
        UG(7,[6 8])=1;
        UG(8,[7 9])=1;
        UG(9,[3 8])=1;
        UG(10,[1 11])=1;
        UG(11,[10 12])=1;
        UG(12,[11 13])=1;
        UG(13,[12 14])=1;
        UG(14,[13 15])=1;
        UG(15,[14 16])=1;
        UG(16,[2 15])=1;
        UG=(UG+UG')>0;
    elseif strcmp(type,'clique')
        % disconnected cliques network
        UG=zeros(16);
        UG(1,2:5)=1;
        UG(2,3:5)=1;
        UG(3,4:5)=1;
        UG(4,5)=1;
        UG(6,7:9)=1;
        UG(7,8:9)=1;
        UG(8,9)=1;
        UG(10,11:12)=1;
        UG(11,12)=1;
        UG(13,14)=1;
        UG=(UG+UG')>0;
    else
        UG=[];
        disp('UNKNOWN TYPE');
    end
end


