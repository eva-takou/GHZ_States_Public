function OUT=Create_GHZ_Sequential(GHZsize)
%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 27, 2023
%--------------------------------------------------------------------------
%
%Script to find optimal entangling gates for creating GHZ states via the sequential
%scheme, based on arXiv:2302.05580. 
%--------------------------------------------------------------------------
%
%Input: GHZsize: size of GHZ state
%Output: OUT: A struct with various fields that stores the several optimal
%             cases together with information about the total evolution.
%--------------------------------------------------------------------------

%======== Parameters ======================================================

[s0,s1,wL] = load_fixed_params;  % Electron's spin projections, and nuclear larmor freq
[A,B]      = Get_HF_Delft;       % HF parameters of 27 spins
Nnuc       = length(A);          % Number of nuclei

[Target_One_Tangle_Tol,Unwanted_One_Tangle_Tol,Infid_Tol,Time_Max,...
 Individual_Time_Max,MaxRes]=Tolerances_Sequential(GHZsize);

%=========================================================================
numSpin           = 1;
Number_Of_Maxima  = 15;
Nuclear_Combs     = nchoosek(1:23,GHZsize-1); %Last 4 spins have B=0 (undetermined by experiment)
%=========================================================================
  
parfor Spin_Indx = 1:23 
    
    times_and_iters{Spin_Indx}=Get_Iters_All_Res(Spin_Indx,MaxRes,wL,A,B,...
                                               s0,s1,Number_Of_Maxima,Individual_Time_Max,...
                                               Unwanted_One_Tangle_Tol) ; 
    
end

%--Test all possible nuclear combinations to select GHZ-1 out of 23 spins--

parfor SpinCombs = 1:size(Nuclear_Combs,1) 
    
    SpinsToTest = Nuclear_Combs(SpinCombs,:);
    [tt,NN]     = Combine_Iters_and_Times(times_and_iters,SpinsToTest,Time_Max);
  
    Times{SpinCombs} = tt; %Times/Iters arranged in terms of spin combinations
    Iters{SpinCombs} = NN; %Times/Iters could contain nan values, for combinations that do not satisfy tolerances
    
end
clear times_and_iters %Free up memory


%------- Remove the nan values --------------------------------------------
to_remove=[];
for SpinCombs=1:size(Nuclear_Combs,1)  
    
   if isnan(Iters{SpinCombs})
       
       to_remove=[to_remove,SpinCombs];
       
   end
    
end

Iters(to_remove)=[];
Times(to_remove)=[];
Nuclear_Combs(to_remove,:)=[];

TIMES = parallel.pool.Constant(Times);       ITERS = parallel.pool.Constant(Iters);
A_HF  = parallel.pool.Constant(A);           B_HF  = parallel.pool.Constant(B);

%=== Find the max # of cases out of all nuclear spin combinations  ========
%   (This is needed to be able to define Max_Case to loop for parfor later)

this_max = size(Iters{1},1);     
LL       = size(Nuclear_Combs,1); %Redefine possible nuclear spin combinations

for ii=2:LL
    
    [new_max,~]=size(Iters{ii});
   
    if new_max>this_max
    
        this_max=new_max; %keep comparing so now Max_Case is this_max
   
    end
    
end

Max_Case = this_max;

%=Get target onetangles for each spin combination & all cases of this comb=

parfor SpinCombs = 1:LL
    
    SpinsToTest       = Nuclear_Combs(SpinCombs,:);   %Target nuclei of this combination 
    [true_max_case,~] = size(Iters{SpinCombs});
    
    for Case=1:Max_Case
        
        if Case<=true_max_case
            
            jj=1;
            
            while jj<=length(SpinsToTest) 
            
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


%==Choose optimal case based on the M-tangle being the cost function=======

CNT=0;

for SpinCombs = 1:LL
    
    case_cnt = 0; 
    CaseNum  = []; 
    Mtangle  = []; %Redefine per spin combination 
    
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
    
    if flag==0 %Found >=1 case within time constraints
        
        CNT       = CNT+1;
        [~,indx]  = max(Mtangle);  %Cost function the M-tangle -- maximal out of all possible cases for this spin combination
        
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
    
    if isempty(ep_unw) && flag==1 %Satisfied one-tangle tolerances
        
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






