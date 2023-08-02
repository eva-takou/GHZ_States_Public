function [byparts]=all_byparts(n)
   %Outputs: 2^n bypartitions in a cell array that need to be taken
   %columnwise from each cell.
    b=1:n;

    byparts = cell(1,n+1); %but 2^n in total when selected columnwise
    for ii=1:n

       byparts{ii}=nchoosek(b,ii).'; 

    end
    %and add+1 which will be the empty set

    byparts{n+1}=[];           
    byparts=circshift(byparts,1); %bring the empty set in the beginning



end
