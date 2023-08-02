function OUT=Create_GHZ_Sequential(GHZsize)


%======== Parameters ======================================================

[s0,s1,wL] = load_fixed_params;  % Electron's spin projections, and nuclear larmor freq
[A,B]      = Get_HF_Delft;       % HF parameters of 27 spins
Nnuc       = length(A);          % Number of nuclei

[Target_One_Tangle_Tol,Unwanted_One_Tangle_Tol,Infid_Tol,Time_Max,Individual_Time_Max,MaxRes]=...
                     Tolerances_Sequential(GHZsize);

%=========================================================================
numSpin           = 1;
Number_Of_Maxima  = 15;
Nuclear_Combs     = nchoosek(1:23,GHZsize-1); %Last 4 spins have B=0 (undetermined by experiment)
%=========================================================================

tic   
parfor Spin_Indx = 1:23 
    
    times_and_iters{Spin_Indx}=Get_Iters_All_Res(Spin_Indx,MaxRes,wL,A,B,...
                                               s0,s1,Number_Of_Maxima,Individual_Time_Max,...
                                               Unwanted_One_Tangle_Tol) ; 
    
end
toc

tic
parfor SpinCombs = 1:size(Nuclear_Combs,1) 
    
    SpinsToTest = Nuclear_Combs(SpinCombs,:);
    [tt,NN]     = Combine_Iters_and_Times(times_and_iters,SpinsToTest,Time_Max);
  
    Times{SpinCombs} = tt;
    Iters{SpinCombs} = NN; 
    
end
toc
clear times_and_iters %Free up memory

%cnt=0;
to_remove=[];
for SpinCombs=1:size(Nuclear_Combs,1)  %Remove any nan values from iterations, times, & nuclear spin combinations
    
   if isnan(Iters{SpinCombs})
       
       to_remove=[to_remove,SpinCombs];
%        cnt                      = cnt+1;
%        Iters_new{cnt}           = Iters{SpinCombs};
%        Times_new{cnt}           = Times{SpinCombs};
%        Nuclear_Combs_new(cnt,:) = Nuclear_Combs(SpinCombs,:);
       
   end
    
end

Iters(to_remove)=[];
Times(to_remove)=[];
Nuclear_Combs(to_remove,:)=[];

TIMES = parallel.pool.Constant(Times);       ITERS = parallel.pool.Constant(Iters);
A_HF  = parallel.pool.Constant(A);           B_HF  = parallel.pool.Constant(B);

%=== Find the max # of cases out of all nuclear spin combinations  ========

this_max = size(Iters{1},1);
LL       = size(Nuclear_Combs,1);

for ii=2:LL
    
    [new_max,~]=size(Iters{ii});
   
    if new_max>this_max
    
        this_max=new_max; %keep comparing so now Max_Case is this_max
   
    end
    
end

Max_Case = this_max;

%====== Get the one-tangles for each case  ================================
tic
parfor SpinCombs = 1:LL
    
    SpinsToTest       = Nuclear_Combs(SpinCombs,:);    
    [true_max_case,~] = size(Iters{SpinCombs});
    
    for Case=1:Max_Case
        
        if Case<=true_max_case
            
            jj=1;
            
            while jj<length(SpinsToTest)
            
                Spin_Indx = SpinsToTest(jj);
                NucSpin   = SubClass_U4Operations(wL,A_HF.Value(Spin_Indx),B_HF.Value(Spin_Indx),s0,s1,numSpin,1,1);   %Res and N doesnt matter since
                NucSpin   = NucSpin.CPMG(TIMES.Value{SpinCombs}(Case,:),ITERS.Value{SpinCombs}(Case,:));
                NucSpin   = NucSpin.Makhlin_Inv; 
                Ep{SpinCombs}{Case}(jj) = NucSpin.Ep/(2/9);
                 
                disp(['Runing case ',num2str(Case),' for Spins ',num2str(SpinsToTest)])
                jj=jj+1;
    
            end
            
        end
        
    end
    
end
toc

%==Choose optimal case based on the M-tangle being the cost function=======

CNT=0;

for SpinCombs = 1:LL
    
    case_cnt = 0; 
    CaseNum  = []; 
    Mtangle  = []; %temp=[];
    
    flag         = 0;
    SpinsToTest  = Nuclear_Combs(SpinCombs,:);
    [Max_Case,~] = size(Iters{SpinCombs});
    
    for Case=1:Max_Case
        
        gate_time = sum(TIMES.Value{SpinCombs}(Case,:).*ITERS.Value{SpinCombs}(Case,:));
        
        if gate_time<=Time_Max %Within time constraints for this spin combination.
            
            case_cnt          = case_cnt+1;
            Mtangle(case_cnt) = prod(Ep{SpinCombs}{Case});
            CaseNum(case_cnt) = Case; 
            
        else
            
            if case_cnt == 0 && Case == Max_Case
                flag=1;
                disp(['For spins ',num2str(SpinsToTest),' could not satisfy time constraints.'])
            end
            
        end
        
    end
    
    if flag==0 %we found at least one case within time cosntraints
        
        CNT       = CNT+1;
        [~,indx]  = max(Mtangle);  %Cost function the M-tangle -- might also consider shortest gate time...
        
        Best_Times{CNT}           = TIMES.Value{SpinCombs}(CaseNum(indx),:);
        Best_Iters{CNT}           = ITERS.Value{SpinCombs}(CaseNum(indx),:);
        Final_NuclearCombs(CNT,:) = SpinsToTest;
        
    end
    
end

[LL,~]=size(Final_NuclearCombs);
%============Completed pre-selection step==================================

%======= Now, check whether the results make sense ========================
disp('-------Entering calculation of tangles in composite evolution.-------')
cnt=0;
for SpinCombs = 1 : LL
    
    times        = Best_Times{SpinCombs};
    iters        = Best_Iters{SpinCombs};
    gate_time    = sum(times.*iters);
    SpinsToTest  = Final_NuclearCombs(SpinCombs,:);
    
    Ep    = zeros(1,Nnuc);  
    n0n1  = zeros(1,Nnuc); 
    n0    = zeros(Nnuc,3);  
    n1    = zeros(Nnuc,3);
    phi0N = zeros(1,Nnuc); 
    phi1N = zeros(1,Nnuc);
    
    for Spin_Indx = 1:Nnuc %evolve all nuclear spins under the composite gate
       
        NucSpin = SubClass_U4Operations(wL,A_HF.Value(Spin_Indx),B_HF.Value(Spin_Indx),s0,s1,numSpin,1,1);   %Res and N doesnt matter 
        NucSpin = NucSpin.CPMG(times,iters);
        NucSpin = NucSpin.Makhlin_Inv;
        NucSpin = NucSpin.Rot_Angles_And_Axes;
     
        Ep(Spin_Indx) = NucSpin.Ep/(2/9);
        
        phi0N(Spin_Indx) = NucSpin.angles{1};
        phi1N(Spin_Indx) = NucSpin.angles{2};
     
        n0(Spin_Indx,:)  = NucSpin.axes{1};
        n1(Spin_Indx,:)  = NucSpin.axes{2};
        n0n1(Spin_Indx)  = dot(NucSpin.axes{1},NucSpin.axes{2});
     
    end
    
    ep_target           = Ep(SpinsToTest);
    ep_unw              = Ep;
    ep_unw(SpinsToTest) = nan;
    ep_unw              = ep_unw(ep_unw>Unwanted_One_Tangle_Tol);
    
    if all(ep_target>Target_One_Tangle_Tol) %For all target nuclei. 
        
        flag=1;
    else
        flag=0;
    end
    
    Unwanted_Nuc = setxor(1:27,SpinsToTest);
    
    if isempty(ep_unw) && flag==1
        
       cnt=cnt+1;

       OUT.Target_Nuclei{cnt}            = SpinsToTest;
       OUT.Unwanted_Nuclei_Names{cnt}    = Unwanted_Nuc;
    
       OUT.A_Target{cnt}                 = A(SpinsToTest);
       OUT.B_Target{cnt}                 = B(SpinsToTest);
    
       OUT.A_Unwanted{cnt}               = A(Unwanted_Nuc);
       OUT.B_Unwanted{cnt}               = B(Unwanted_Nuc);
    
       OUT.Opt_Unit_Times{cnt}             = times;
       OUT.Opt_Iters{cnt}                  = iters;
    
       OUT.Total_Time{cnt}                 = gate_time;
       OUT.EP_Target{cnt}                  = Ep(SpinsToTest);
       OUT.EP_Unwanted{cnt}                = Ep(Unwanted_Nuc);
    
       OUT.phi0N{cnt}                      = phi0N;
       OUT.phi1N{cnt}                      = phi1N;
       OUT.n0{cnt}                         = n0;
       OUT.n1{cnt}                         = n1;
       OUT.n0n1{cnt}                       = n0n1;
       
    else
        
        disp('Excluded additional cases based on tolerances of one-tangles.')
       
    end
    
end

%======= Finally, get the gate errors. ====================================

Max_Cases = cnt;

disp(['-----------GHZ instances:',num2str(Max_Cases),'-----------'])
disp('-----------Entering calculations of gate errors.-----------')

OUT_parfor = parallel.pool.Constant(OUT);

parfor Cases  = 1 : Max_Cases
    
    Unwanted_Nuclei = OUT_parfor.Value.Unwanted_Nuclei_Names{Cases};
    ep_unw          = OUT_parfor.Value.EP_Unwanted{Cases}; %has same length as Unwanted_Nuclei
    
    Unwanted_Nuclei(ep_unw<1e-8) = nan;
    Unwanted_Nuclei              = Unwanted_Nuclei(~isnan(Unwanted_Nuclei));
    
    phi0_Unw        = OUT_parfor.Value.phi0N{Cases}(Unwanted_Nuclei);
    phi1_Unw        = OUT_parfor.Value.phi1N{Cases}(Unwanted_Nuclei);
    n0_Unw          = OUT_parfor.Value.n0{Cases}(Unwanted_Nuclei,:);
    n1_Unw          = OUT_parfor.Value.n1{Cases}(Unwanted_Nuclei,:);
    temp            = SubClass_Ent_and_Fid;
    Infid(Cases)    = temp.Gate_Infid(length(Unwanted_Nuclei)+GHZsize-1,...
                      GHZsize-1,phi0_Unw,phi1_Unw,n0_Unw,n1_Unw);
    
end
clear OUT_parfor

%======= Post-select based on gate error ==================================

to_remove=[];
for Cases = 1:Max_Cases
    
   OUT.Infid{Cases} =  Infid(Cases);
    
   if Infid(Cases).Infid>Infid_Tol
       
       to_remove=[to_remove,Cases];
       
   end
   
end


fields = fieldnames(OUT);

        
for jj=1:length(fields)
    
   temp = OUT.(fields{jj});
   temp(to_remove)=[];
   OUT.(fields{jj})=temp;
    
end




end

function [t_opt,N]=Get_Iters(Spin_Indx,Res,wL,A,B,s0,s1,Number_Of_Maxima,Time_Max,Unwanted_One_Tangle_Tol)
%Output all possible Iters and optimal time of a particular resonance and
%nuclear spin such that the unwanted one-tangles are below the tolerance.
%Input: Spin_Indx \in [1,27]
%       Res: resonance indx
%       wL,A,B: Larmor frequency and HF parameters of 27 nuclei
%       s0,s1: Electron's projections
%       Number_Of_Maxima: how many maxima of one-tangle to keep 
%       Time_Max: Time constraint
%       Unwanted_One_Tangle_Tol: Tolerance of unwanted one-tangles.
%Output: optimal unit time and iterations N.
Nnuc    = 27;
TIME    = Resonance_Time(Res,Spin_Indx,wL,A,B,s0,s1);
tt      = TIME-0.2:5e-3:TIME+0.2;
numSpin = 1;
n0n1    = zeros(1,length(tt));

for ii=1:length(tt) %Find when the dot product==-1
    
NucSpin = SubClass_U4Operations(wL,A(Spin_Indx),B(Spin_Indx),s0,s1,numSpin,Res,1);  %1 iteration
NucSpin  = NucSpin.CPMG(tt(ii));
NucSpin  = NucSpin.Rot_Angles_And_Axes;
n0n1(ii) = dot(NucSpin.axes{1},NucSpin.axes{2} );

end

[~,indx]=min(n0n1);    t_opt = tt(indx);
clear tt
clear n0n1
NucSpin = SubClass_U4Operations(wL,A(Spin_Indx),B(Spin_Indx),s0,s1,numSpin,Res,1);  %1 iteration)
NucSpin = NucSpin.CPMG(t_opt);
NucSpin = NucSpin.Expected_Maxima(Number_Of_Maxima);
N       = NucSpin.Nmax;

for ii=1:length(N)
    
   if ~isreal(N(ii)) || N(ii)<0 
   
       error('Error in the evaluation of iterations that maximize one-tangles.')
       
   end
    
end

N = N(N*t_opt<Time_Max); %Restrict N such that we are within time limits.

if isempty(N) 
    
    t_opt = nan;
    N     = nan;
    
    return
    
end

%Remove N where epunwanted>unwanted one-tangle tolerance 
for ii=1:length(N)
    
    spin=1;
    
    while spin<=Nnuc
         
        NucSpin = SubClass_U4Operations(wL,A(spin),B(spin),s0,s1,numSpin,Res,N(ii));  
        NucSpin = NucSpin.CPMG(t_opt);
        NucSpin = NucSpin.Makhlin_Inv;
        Ep      = NucSpin.Ep/(2/9);
        
        if spin~=Spin_Indx && Ep>Unwanted_One_Tangle_Tol
            
            N(ii)=nan;
            break
            
        end
        spin=spin+1;
    end
    
end

N = N(~isnan(N));

if isempty(N)
    
    t_opt = nan;
    N     = nan;
end

                   
end


function OUT=Get_Iters_All_Res(Spin_Indx,Res_Max,wL,A,B,s0,s1,Number_Of_Maxima,Time_Max,Unwanted_One_Tangle_Tol)
%Output all possible Iters and optimal times for up to a max resonance and
%for a specific nuclear spin labeled by Spin_Indx.
%Input: Spin_Indx \in [1,27]
%       Res_Max: Maximum resonance to search for
%       wL,A,B: Larmor frequency and HF parameters of 27 nuclei
%       s0,s1: Electron's projections
%       Number_Of_Maxima: how many maxima of one-tangle to keep 
%       Time_Max: Time constraint
%       Unwanted_One_Tangle_Tol: Tolerance of unwanted one-tangles.
%Output: Optimal unit times and iterations.
cnt=0;
CNT=0;

for Res=1:Res_Max
    
   [t,N]=Get_Iters(Spin_Indx,Res,wL,A,B,s0,s1,Number_Of_Maxima,Time_Max,Unwanted_One_Tangle_Tol) ;
    
    if ~isnan(t)
        
       for ll=1:length(N)
       
           cnt            = cnt+1; 
           OUT.Times(cnt) = t;
           OUT.Iters(cnt) = N(ll);
       
       end
       
    else
        
        CNT=CNT+1;
       
    end
                   
    
end

if CNT==Res_Max %We didnt satisfy the various constraints. No acceptable cases.
    
   OUT=nan; 
end


end


function [Times,Iters]=Combine_Iters_and_Times(input_cell,Spin_Indices,Time_Max)
%Use this to get all combinations of times and iterations.

for ii=1:length(Spin_Indices)
    
    if isstruct(input_cell{Spin_Indices(ii)})
    
        times_to_comb{ii} = input_cell{Spin_Indices(ii)}.Times;
        iters_to_comb{ii} = input_cell{Spin_Indices(ii)}.Iters;
        
    else %not a good spin combination
        
        Times=nan;
        Iters=nan;
        return
    end
    
end

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
           
           Times(cnt,ii) = combinations_T{ii}(kk); 
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

