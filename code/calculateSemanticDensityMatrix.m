function [ crosscount ] = calculateSemanticDensityMatrix( MatrixProfile, MPindex, k ,m)
%UNTITLED Summary of this function goes here
%   This function extracts k-chaine of simillar patterns for each
%   subsequence. input:
%   MatrixProfile, MPindex
%   Outputs:
%   extracted segment indices
dontcare=length(MatrixProfile)+1;
%threshold = mean(MatrixProfile,'omitnan')+3*std(MatrixProfile,'omitnan');
threshold = 2*max(MatrixProfile);
if (size(MatrixProfile,1)==1)
    MatrixProfile = MatrixProfile';
end

if (size(MPindex,1)==1)
    MPindex = MPindex';
end

ArcSet=num2cell(MPindex);
ArcCost=num2cell(MatrixProfile);

lastArcSet=MPindex;
lastArcCost=MatrixProfile;
for i=1:k
    [ArcSet,ArcCost,lastArcSet,lastArcCost]=ExtractNewArcSet...
        (MatrixProfile,MPindex,ArcSet,ArcCost,threshold,lastArcSet,lastArcCost,m);
    if (sum(lastArcSet(:)<dontcare)==0)
        break;
    end
    sum(lastArcSet(:)<dontcare);
end
nnmark=zeros(1,length(MatrixProfile));
%MAX AND MIN
MIN=cellfun(@(x) min(x),ArcSet);
MAX=cellfun(@(x) max(x),ArcSet);
totmin=min(abs(MIN-[1:size(MIN)]'));
totmax=max(abs(MAX-[1:size(MAX)]'));

for j=1:size(ArcSet,1)
    list = ArcSet{j};
    list_len = length(list);
    % find the number of arcs crossing over each index i
    %nnmark2=zeros(1,profile_len);
    for i=1:list_len
        small=min(j,list(i));
        large=max(j,list(i));
        len=large-small;
        %nnmark(small:large)=nnmark(small:large)+1; %SHOHREH
        nnmark(small:large)=nnmark(small:large)+(1-((len-totmin)/(totmax-totmin))); %SHOHREH
    end
end
crosscount = nnmark;
end

function [ArcSet,ArcCost,lastArcSet,lastArcCost]=ExtractNewArcSet...
    (MatrixProfile,MPindex,ArcSet,ArcCost,threshold,lastArcSet,lastArcCost,m)

dontcare=length(MatrixProfile)+1;
InitialArcs = lastArcSet;
InitialArcs(lastArcCost(:)>threshold)=dontcare;

temp = [InitialArcs; dontcare];

newArcs = temp(InitialArcs);
tempArcCost(1:dontcare-1) =threshold+1;
tempArcSet(1:dontcare-1) =dontcare;
for i = 1 : size(ArcSet,1)
    if(newArcs(i)~=dontcare && (newArcs(i)>i+m/4 || newArcs(i)<i-m/4))
            ArcSet{i}=[ArcSet{i},newArcs(i)];
            ArcCost{i} = [ArcCost{i},lastArcCost(i)+ lastArcCost(lastArcSet(i))];
            tempArcCost(i)=lastArcCost(i)+ lastArcCost(lastArcSet(i));
            tempArcSet(i)=newArcs(i);
        
    else
        tempArcCost(i) = threshold+1;
        tempArcSet(i) = dontcare;
    end
end
lastArcCost=[tempArcCost]';
lastArcSet=[tempArcSet]';
end