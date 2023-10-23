function [V,d]=sorted_evals_evecs(H)
%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
%
%Sort the evals of H, and normalize the e-vecs of H.
%--------------------------------------------------------------------------

[dim1,~]=size(H);

BB=eye(dim1);

[V,D] = eig(H,BB,'qz');  

d      = diag(sort(real(diag(D)),'descend'));
d      = diag(d);

[~,ind]=sort(real(diag(D)),'descend'); %store the indices of which columns the sorted evals come from

V=V(:,ind); %arange the evecs according to order of evals

[~,col]=size(V);

for ii=1:col
   V(:,ii)=V(:,ii)/norm(V(:,ii)); 
    
end


end
