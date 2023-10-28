function OUT=Get_ep_for_fixed_k(resonance,wL,s0,s1,A,B,Time_Max,...
                       Unwanted_One_Tangle_Tol,Target_One_Tangle_Tol,GHZsize,Length_Of_Time_Vals)

%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 28, 2023
%--------------------------------------------------------------------------
%           
%Function to get the one-tangles, times and iterations where GHZ-1 
%one-tangles become maximal (at the same time and iterations).
%
%Input: resonance: resonance index
%       wL: Larmor frequency of nuclei
%       A/B: Parallel/Perpendicular HF components
%       Time_Max: Max time of sequence
%       Unwanted_One_Tangle_Tol/Target_One_Tangle_Tol: Tolerances for
%       values one-tangles
%       GHZsize: size of GHZ state to create
%       Length_of_Time_Vals: # of total values of time to sample over (in a \pm 0.25 resonance range)
%
%Output: A struct with possible times, iterations, one-tangles, and indices of target
%        nuclei.


Nnuc       = length(A);
unit_times = zeros(1,Nnuc);
tt         = cell(1,24);

for fixed_Nuc = 1:24

    unit_times(fixed_Nuc) = Resonance_Time(resonance,fixed_Nuc,wL,A,B,s0,s1);
    tmin                  = (1-0.25)    * (unit_times(fixed_Nuc)); 
    tmax                  = (1+0.25)    * (unit_times(fixed_Nuc)); 
    tt{fixed_Nuc}         = linspace(tmin,tmax,Length_Of_Time_Vals);
    
end

tt = cell2mat(tt);
tt = tt(:);
Lt = length(tt);

NMAX = round(Time_Max/min(tt));

Times_Iters = allcomb(tt,1:NMAX);
Times       = Times_Iters(:,1);
Iters       = round(Times_Iters(:,2));

clear Times_Iters;

disp(['Times I will loop for:',num2str(Lt),'...'])

%- Evolve all  nuclei for each time and collect rot angles and dot prod---

parfor ii=1:Lt
    
    for Spin_Indx=1:Nnuc
        
       numSpin            = 1;    
       NucSpin            = SubClass_U4Operations(wL,A(Spin_Indx),B(Spin_Indx),s0,s1,numSpin,resonance,1);    
       NucSpin            = NucSpin.CPMG(tt(ii));      
       NucSpin            = NucSpin.Rot_Angles_And_Axes;
       phi0(ii,Spin_Indx) = NucSpin.angles{1};
       phi1(ii,Spin_Indx) = NucSpin.angles{2};
       n0n1(ii,Spin_Indx) = dot(NucSpin.axes{1},NucSpin.axes{2});
        
    end
    
end
disp('Done.')

disp('Collecting one-tangles...')
%-Evaluate one tangle as a function of N ----------------------------------
parfor ii=1:Lt
    
    for NN=1:NMAX
     
        phi0N        = acos(cos(NN.*phi0(ii,:)));
        phi1N        = acos(cos(NN.*phi1(ii,:)));
        G1          = ( cos(phi0N/2).*cos(phi1N/2) + n0n1(ii,:).*sin(phi0N/2).*sin(phi1N/2) ).^2;
        Ep(NN,ii,:) = 1-G1; 
     
    end
end
disp('Done.')

%------- This is test for alternative evaluation of Ep --------------------
% parfor ii=1:Lt
%     
%     for NN=1:4
%        
%         for jj=1:Nnuc
%             
%             NucSpin  = SubClass_U4Operations(wL,A(jj),B(jj),s0,s1,1,1,NN);    
%             NucSpin  = NucSpin.CPMG(tt(ii));
%             NucSpin  = NucSpin.Makhlin_Inv;
%             
%             Ep_Alt(NN,ii,jj)=NucSpin.Ep/(2/9);
%             
%         end
%         
%     end
%     
% end
%--------------------------------------------------------------------------

clear phi0; clear phi1; clear n0n1;

Ep = reshape(Ep,[length(Iters),Nnuc]); % [length(Iters) x length(Times)]x Nnuc
                                       
%for some N we exceeded the gate time
indices =  Times.*Iters<Time_Max;
Times   =  Times(indices);
Iters   =  Iters(indices);
Ep      =  Ep(indices,:);

Max_Case  = length(Times);
cnt       = 0;

%--------- Another test to make sure that reshape worked correctly --------

% parfor jj=1:length(Times)
%     
%     for k=1:Nnuc
%         
%         test=SubClass_U4Operations(wL,A(k),B(k),s0,s1,1,1,Iters(jj));
%         test=test.CPMG(Times(jj));
%         test=test.Makhlin_Inv;
%         
%         ep_test=test.Ep/(2/9);
%         
%         if abs(ep_test-Ep(jj,k))>1e-8
%             
%            error('Incorrect reshaping') 
%             
%         end
%         
%     end
% end
%--------------------------------------------------------------------------


disp('Checking acceptable cases based on one-tangles tolerances...')
parfor Case = 1 : Max_Case
   
    ep_All    = Ep(Case,:);
    ep_T      = ep_All(ep_All>Target_One_Tangle_Tol);
    ep_U      = ep_All(ep_All<Target_One_Tangle_Tol);
 
    if length(ep_T)==GHZsize-1 && isempty(ep_U(ep_U>Unwanted_One_Tangle_Tol)) 
        %If this is not satisfied, then the parfor will put empty thing on OUT_new2{Case}                       
        
        cnt=cnt+1;
        OUT_new2{Case}.Target_Nuclei = find(ep_All>Target_One_Tangle_Tol);
        OUT_new2{Case}.Ep            = Ep(Case,:)
        OUT_new2{Case}.Times         = Times(Case);
        OUT_new2{Case}.Iters         = Iters(Case); 
    end
    
end


disp('Done.')

if cnt==0   %we failed
   
   disp(['Found no acceptable case for k=',num2str(resonance)]) 
   
   OUT=nan; 
   return 
   
else
    clear Iters; clear Times; clear Ep;
    disp(['Found some acceptable cases for k=',num2str(resonance)])
   
end

cnt=0;
for Case=1:length(OUT_new2)
    
    if isstruct(OUT_new2{Case})
        cnt=cnt+1;
        OUT.Ep(cnt,:)            = OUT_new2{Case}.Ep;
        OUT.Times(cnt)           = OUT_new2{Case}.Times;
        OUT.Iters(cnt)           = OUT_new2{Case}.Iters;
        OUT.Target_Nuclei(cnt,:) = OUT_new2{Case}.Target_Nuclei;
    end
    
    
end


disp('-----------------------------------------------------')


end

