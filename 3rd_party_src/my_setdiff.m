function C = my_setdiff(A,B)
%
%Nick (2023). setdiff (https://www.mathworks.com/matlabcentral/fileexchange/23172-setdiff), MATLAB Central File Exchange. Retrieved October 23, 2023.
%
% MYSETDIFF Set difference of two sets of positive integers (much faster than built-in setdiff)
% C = my_setdiff(A,B)
% C = A \ B = { things in A that are not in B }
% 
%if isempty(A)
 %   C = [];
  %  return;
if isempty(B)
    C = A;
    return; 
else % both non-empty
    bits = zeros(1, max(max(A), max(B)));
    bits(A) = 1;
    bits(B) = 0;
    C = A(logical(bits(A)));
end         




end
