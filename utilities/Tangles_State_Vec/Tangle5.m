function tau=Tangle5(psi)


temp_1=0; temp_2=0; temp_3=0; temp_4=0;  temp_5=0;

n=5;
bitCombs = dec2bin(0:2^n-1) - '0'; %each row is a combination of bitstrings


for ii=1:length(bitCombs)
    
    for jj=1:length(bitCombs)
        for kk=1:length(bitCombs)
            for mm=1:length(bitCombs)
   
               fixed=to_ket(bitCombs(ii,:),psi)*...
                     to_ket(bitCombs(jj,:),psi)*...
                     to_ket(bitCombs(kk,:),psi)*...
                     to_ket(bitCombs(mm,:),psi);
  

               
             fixed1=  ee(bitCombs(ii,2),bitCombs(jj,2))*ee(bitCombs(ii,3),bitCombs(jj,3))*ee(bitCombs(ii,4),bitCombs(jj,4))*ee(bitCombs(ii,5),bitCombs(jj,5))*...
                      ee(bitCombs(kk,2),bitCombs(mm,2))*ee(bitCombs(kk,3),bitCombs(mm,3))*ee(bitCombs(kk,4),bitCombs(mm,4))*ee(bitCombs(kk,5),bitCombs(mm,5))*...
                      fixed;

             fixed2=  ee(bitCombs(ii,1),bitCombs(jj,1))*ee(bitCombs(ii,3),bitCombs(jj,3))*ee(bitCombs(ii,4),bitCombs(jj,4))*ee(bitCombs(ii,5),bitCombs(jj,5))*...
                      ee(bitCombs(kk,1),bitCombs(mm,1))*ee(bitCombs(kk,3),bitCombs(mm,3))*ee(bitCombs(kk,4),bitCombs(mm,4))*ee(bitCombs(kk,5),bitCombs(mm,5))*...
                      fixed;

             fixed3=  ee(bitCombs(ii,1),bitCombs(jj,1))*ee(bitCombs(ii,2),bitCombs(jj,2))*ee(bitCombs(ii,4),bitCombs(jj,4))*ee(bitCombs(ii,5),bitCombs(jj,5))*...
                      ee(bitCombs(kk,1),bitCombs(mm,1))*ee(bitCombs(kk,2),bitCombs(mm,2))*ee(bitCombs(kk,4),bitCombs(mm,4))*ee(bitCombs(kk,5),bitCombs(mm,5))*...
                      fixed;

              fixed4= ee(bitCombs(ii,1),bitCombs(jj,1))*ee(bitCombs(ii,2),bitCombs(jj,2))*ee(bitCombs(ii,3),bitCombs(jj,3))*ee(bitCombs(ii,5),bitCombs(jj,5))*...
                      ee(bitCombs(kk,1),bitCombs(mm,1))*ee(bitCombs(kk,2),bitCombs(mm,2))*ee(bitCombs(kk,3),bitCombs(mm,3))*ee(bitCombs(kk,5),bitCombs(mm,5))*...
                      fixed;

              fixed5= ee(bitCombs(ii,1),bitCombs(jj,1))*ee(bitCombs(ii,2),bitCombs(jj,2))*ee(bitCombs(ii,3),bitCombs(jj,3))*ee(bitCombs(ii,4),bitCombs(jj,4))*...
                      ee(bitCombs(kk,1),bitCombs(mm,1))*ee(bitCombs(kk,2),bitCombs(mm,2))*ee(bitCombs(kk,3),bitCombs(mm,3))*ee(bitCombs(kk,4),bitCombs(mm,4))*...
                      fixed;
               
               
 temp_1=temp_1+fixed1*ee(bitCombs(ii,1),bitCombs(kk,1))*ee(bitCombs(jj,1),bitCombs(mm,1));
 temp_2=temp_2+fixed2*ee(bitCombs(ii,2),bitCombs(kk,2))*ee(bitCombs(jj,2),bitCombs(mm,2));
 temp_3=temp_3+fixed3*ee(bitCombs(ii,3),bitCombs(kk,3))*ee(bitCombs(jj,3),bitCombs(mm,3));
 temp_4=temp_4+fixed4*ee(bitCombs(ii,4),bitCombs(kk,4))*ee(bitCombs(jj,4),bitCombs(mm,4));
 temp_5=temp_5+fixed5*ee(bitCombs(ii,5),bitCombs(kk,5))*ee(bitCombs(jj,5),bitCombs(mm,5));
             
               
               
            end
        end
    end
    
    
    
    
end

 
                                                            

                   

temp_1=abs(temp_1);
temp_2=abs(temp_2);
temp_3=abs(temp_3);
temp_4=abs(temp_4);
temp_5=abs(temp_5);


tau = 2/5*(temp_1+temp_2+temp_3+temp_4+temp_5);





end


function out=to_ket(indices,a)




indices=flip(indices);
num=0;
for ii=0:length(indices)-1
   
    indx = indices(ii+1);
    
    num = num+ indx*2^(ii) ;
    
    
end

%the input is the state |ijklm>

%num = m*2^0 + l*2^1 + k*2^2 + j*2^3 + i*2^4

%num = num + 1;

out = a(num+1);








end


function out=ee(l,m)


if l==m %0,0 or 1,1
    out=0;
elseif l==0 && m==1
    out=1;
elseif l==1 && m==0
    out=-1;
    
end


end



