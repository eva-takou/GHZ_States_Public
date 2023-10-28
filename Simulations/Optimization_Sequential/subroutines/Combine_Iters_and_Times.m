function [Times,Iters]=Combine_Iters_and_Times(input_cell,Spin_Indices,Time_Max)
%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 27, 2023
%--------------------------------------------------------------------------
%
%Use this to get all combinations of times and iterations, which are
%acceptable from one-tangle constraints.
%
%Input: input_cell: A cell array where each entry contains times and
%                   iterations that maximize target one-tangle and minimize unwanted
%                   one-tangle. The input is from Get_Iters_All_Res
%       Spin_Indices: To select particular entries from the cell array
%       (particular spins, each with some realizations of times and
%       iterations)
%       Time_Max: The maximum time of the total composite sequence
%
%Output: For the spin indices we chose, we form the combinations of t_j and
%        N_j to perform the sequential gates.
%        Times/Iters: Each is an array of possible cases (rows) x number of spins to control (columns). 


cell_to_test = input_cell(Spin_Indices); %Get the times/iters for the spins we want to sequential entangle with the electron


% -----Make sure that this spin combination has accepted times/iters-------
flag=false(1,length(cell_to_test));

for jj=1:length(cell_to_test)
    
    if isstruct(cell_to_test{jj})
       
        flag(jj)=true;
        
    end
    
end

if all(flag)
    
    for ii=1:length(cell_to_test)
        
        times_to_comb{ii} = cell_to_test{ii}.Times;
        iters_to_comb{ii} = cell_to_test{ii}.Iters;
        
    end
    
else
    
    Times=nan;
    Iters=nan;
    return
    
    
end
%--------------------------------------------------------------------------

%Use this method to get all possible combinations 

combinations_T = cell(1, numel(times_to_comb)); 
[combinations_T{:}] = ndgrid(times_to_comb{:});
combinations_T = cellfun(@(x) x(:), combinations_T,'uniformoutput',false); 

combinations_N = cell(1, numel(iters_to_comb)); %set up the varargout result
[combinations_N{:}] = ndgrid(iters_to_comb{:});
combinations_N = cellfun(@(x) x(:), combinations_N,'uniformoutput',false);

clear times_to_comb %to free some memory
clear iters_to_comb %to free some memory

%now from each cell of combinations_T and combinations_N take 1 element and
%test the total gate time

max_case = length(combinations_T{1});
cnt=0;

for kk=1:max_case
    
    total_time = 0;
   
   for ii=1:length(combinations_T)
       
       total_time  = total_time + combinations_T{ii}(kk)*combinations_N{ii}(kk);
       
   end
   
   if total_time < Time_Max
       
        cnt=cnt+1;
        
       for ii=1:length(combinations_T)
           
           Times(cnt,ii) = combinations_T{ii}(kk); %rows are realizations, columns correspond to time/iter for control of each spin
           Iters(cnt,ii) = combinations_N{ii}(kk);
           
       end
        
   end
end

if cnt==0 %We didnt find acceptable cases.
    
    Times=nan;
    Iters=nan;
    return
    
else
    
    [max_case,~]=size(Times);
    
    for ii=1:max_case
        disp(['Total time:',num2str(sum(Times(ii,:).*Iters(ii,:)))])
    end    
    
end
      
          
end
