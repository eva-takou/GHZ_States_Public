function S=Create_GHZ_Multispin(GHZsize)
%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 28, 2023
%--------------------------------------------------------------------------
%
%GHZsize: Size of the GHZ state (>=4).
%
%
%Script to find optimal entangling gates for creating GHZ states via the
%multispin scheme, based on arXiv:2302.05580. 
%--------------------------------------------------------------------------
%
%Input: GHZsize: size of GHZ state
%Output: S: A struct with various fields that stores the several optimal
%             cases together with information about the total evolution.
%--------------------------------------------------------------------------


if GHZsize==3
   
    %error('This script is for generation of GHZ_M states with M>3.')

end

%============= Set the basic parameters ====================================
[s0,s1,wL] = load_fixed_params;  % Electron's spin projections, and nuclear larmor freq
[A,B]      = Get_HF_Delft;       % HF parameters of 27 spins
Nnuc       = length(A);          % Number of nuclei

[Target_One_Tangle_Tol,Unwanted_One_Tangle_Tol,Infid_Tol,Time_Max,~,Time_Vals]=Tolerances_Multispin(GHZsize);

%============== Run for multiple resonances (up to 10) ===================

OUT{1}  =Get_ep_for_fixed_k(1,wL,s0,s1,A,B,Time_Max,...
                       Unwanted_One_Tangle_Tol,Target_One_Tangle_Tol,GHZsize,Time_Vals);
disp('Finished with 1st resonance')    
    
OUT{2}  = Get_ep_for_fixed_k(2,wL,s0,s1,A,B,Time_Max,...
                       Unwanted_One_Tangle_Tol,Target_One_Tangle_Tol,GHZsize,Time_Vals);
disp('Finished with 2nd resonance')
OUT{3}  = Get_ep_for_fixed_k(3,wL,s0,s1,A,B,Time_Max,...
                       Unwanted_One_Tangle_Tol,Target_One_Tangle_Tol,GHZsize,Time_Vals);
disp('Finished with 3rd resonance')
OUT{4}  = Get_ep_for_fixed_k(4,wL,s0,s1,A,B,Time_Max,...
                       Unwanted_One_Tangle_Tol,Target_One_Tangle_Tol,GHZsize,Time_Vals);
disp('Finished with 4th resonance')
OUT{5}  = Get_ep_for_fixed_k(5,wL,s0,s1,A,B,Time_Max,...
                       Unwanted_One_Tangle_Tol,Target_One_Tangle_Tol,GHZsize,Time_Vals);
disp('Finished with 5th resonance')
OUT{6}  = Get_ep_for_fixed_k(6,wL,s0,s1,A,B,Time_Max,...
                       Unwanted_One_Tangle_Tol,Target_One_Tangle_Tol,GHZsize,Time_Vals);
disp('Finished with 6th resonance')
OUT{7}  = Get_ep_for_fixed_k(7,wL,s0,s1,A,B,Time_Max,...
                       Unwanted_One_Tangle_Tol,Target_One_Tangle_Tol,GHZsize,Time_Vals);
disp('Finished with 7th resonance')
OUT{8}  = Get_ep_for_fixed_k(8,wL,s0,s1,A,B,Time_Max,...
                       Unwanted_One_Tangle_Tol,Target_One_Tangle_Tol,GHZsize,Time_Vals);
disp('Finished with 8th resonance')
OUT{9}  = Get_ep_for_fixed_k(9,wL,s0,s1,A,B,Time_Max,...
                       Unwanted_One_Tangle_Tol,Target_One_Tangle_Tol,GHZsize,Time_Vals);
disp('Finished with 9th resonance')
OUT{10} = Get_ep_for_fixed_k(10,wL,s0,s1,A,B,Time_Max,...
                       Unwanted_One_Tangle_Tol,Target_One_Tangle_Tol,GHZsize,Time_Vals);
disp('Finished with 10th resonance')

%------- Remove nan cases -------------------------------------------------
cnt=0;
for ii=1:length(OUT)
    
    if isstruct(OUT{ii}) %Not nan-found some acceptable cases.
        cnt=cnt+1;
        TEMP{cnt}=OUT{ii};
    end
    
end

OUT = cell2mat(TEMP);

disp('unfolding...')  

OUT_new.Times         = [];
OUT_new.Iters         = [];
OUT_new.Ep            = [];
OUT_new.Target_Nuclei = [];

for ii=1:cnt
    
    OUT_new.Times         = [OUT_new.Times,OUT(ii).Times];
    OUT_new.Iters         = [OUT_new.Iters,OUT(ii).Iters];
    OUT_new.Ep            = [OUT_new.Ep;OUT(ii).Ep];
    OUT_new.Target_Nuclei = [OUT_new.Target_Nuclei;OUT(ii).Target_Nuclei];
end

disp('Done.')

clear OUT %to free memory

NuclearCombs = nchoosek(1:24,GHZsize-1);
[L,~]        = size(NuclearCombs);
Max_Case     = length(OUT_new.Times);
S            = cell(1,L);

disp('Aranging acceptable cases based on unique nuclear combinations...')
parfor comb=1:L
    
    cnt=0;
    Spins_To_Test = NuclearCombs(comb,:);
    
    for Case = 1:Max_Case
    
        if all(OUT_new.Target_Nuclei(Case,:)==Spins_To_Test) %CHECK ALL HERE
            
            cnt=cnt+1;
            S{comb}.Target_Nuclei      = Spins_To_Test;
            S{comb}.Times(cnt)         = OUT_new.Times(Case);
            S{comb}.Iters(cnt)         = OUT_new.Iters(Case);
            S{comb}.Ep(cnt,:)          = OUT_new.Ep(Case,:);
            
        end
        
    end
    
end
disp('Done.')

clear OUT_new %to free memory

%-- There could be spin combs that are empty so remove them ---------------
cnt=0;
to_remove=[];

for ii=1:length(S)
    
    if isstruct(S{ii})
        cnt=cnt+1;
        
    else
        to_remove=[to_remove,ii];
        
    end
    
end

S(to_remove)=[];

disp('===================================================================')
disp(['Total unique spin combinations:',num2str(cnt)])
disp('===================================================================')

%For each unique combination select the optimal (based on maximal M-tangle).

for comb=1:length(S)
    
    Max_Case  = length(S{comb}.Times );
        
    M_Tangle=zeros(1,Max_Case);
    
    for jj=1:Max_Case

        M_Tangle(jj) = prod(S{comb}.Ep(jj,S{comb}.Target_Nuclei)); %S{comb}.Ep is Case x Nnuc array

    end
    
    [~,indx] = max(M_Tangle,[],'all','linear');

    opt_N = S{comb}.Iters(indx);
    opt_t = S{comb}.Times(indx);
    
    S{comb}.Unwanted_Nuclei = setxor(1:27,S{comb}.Target_Nuclei);

    S{comb}.Times         = opt_t;
    S{comb}.Iters         = opt_N;
    S{comb}.Total_Time    = opt_t*opt_N;

    S{comb}.Ep             = S{comb}.Ep(indx,:);
    S{comb}.Ep_Target      = S{comb}.Ep(S{comb}.Target_Nuclei);
    S{comb}.Ep_Unwanted    = S{comb}.Ep(S{comb}.Unwanted_Nuclei);

    S{comb}.A_Target        = A(S{comb}.Target_Nuclei);
    S{comb}.B_Target        = B(S{comb}.Target_Nuclei);

    S{comb}.A_Unwanted      = A(S{comb}.Unwanted_Nuclei);
    S{comb}.B_Unwanted      = B(S{comb}.Unwanted_Nuclei);
    
end

%-------- Test that in each case we evaluated correctly ------------------
% for jj=1:length(S)
%     
%     for k=1:Nnuc
%    
%         spin=SubClass_U4Operations(wL,A(k),B(k),s0,s1,1,1,S{jj}.Iters);
%         spin=spin.CPMG(S{jj}.Times);
%         spin=spin.Makhlin_Inv;
%         ep_Test(k)=spin.Ep/(2/9);
%     
%     end
%     
%     if any(abs(ep_Test-S{jj}.Ep)>1e-8)
%        
%         error('The entangling power was not evaluated correctly')
%         
%     end
%     
%     if any(abs(ep_Test(S{jj}.Target_Nuclei)-S{jj}.Ep_Target)>1e-8)
%         
%         error('The entangling power was not evaluated correctly')
%     end
%     
%     if any(abs(ep_Test(S{jj}.Unwanted_Nuclei)-S{jj}.Ep_Unwanted)>1e-8)
%         
%         error('The entangling power was not evaluated correctly')
%     end
%     
% end
%--------------------------------------------------------------------------

%=================== Get the evolution of all spins =====================

Final_Cases = length(S);

for ii=1:Final_Cases
    
    phi0N = zeros(1,Nnuc);
    phi1N = zeros(1,Nnuc);
    n0    = zeros(Nnuc,3);
    n1    = zeros(Nnuc,3);
    n0n1  = zeros(1,Nnuc);
    
    for Spin_Indx = 1 : Nnuc
        
        numSpin   = 1;    
        NucSpin = SubClass_U4Operations(wL,A(Spin_Indx),B(Spin_Indx),s0,s1,numSpin,1,S{ii}.Iters);                               
        NucSpin = NucSpin.CPMG(S{ii}.Times);
        NucSpin = NucSpin.Rot_Angles_And_Axes;

        phi0N(Spin_Indx) = NucSpin.angles{1};
        phi1N(Spin_Indx) = NucSpin.angles{2};
        n0(Spin_Indx,:)  = NucSpin.axes{1};
        n1(Spin_Indx,:)  = NucSpin.axes{2};
        n0n1(Spin_Indx)  = dot(NucSpin.axes{1},NucSpin.axes{2});
        
    end
    
        S{ii}.phi0N  = phi0N;
        S{ii}.phi1N  = phi1N;
        S{ii}.n0     = n0;
        S{ii}.n1     = n1;
        S{ii}.n0n1   = n0n1;
        
end
    
%======= Finaly get the gate errors =====================================

disp('-----------Entering calculations of gate errors.-----------')

S_parfor = parallel.pool.Constant(S);

parfor ii=1:Final_Cases
    
    Unwanted_Nuclei = S_parfor.Value{ii}.Unwanted_Nuclei;
    ep_unw          = S_parfor.Value{ii}.Ep_Unwanted; %has same length as Unwanted_Nuclei
    
    Unwanted_Nuclei(ep_unw<1e-8) = nan;
    Unwanted_Nuclei              = Unwanted_Nuclei(~isnan(Unwanted_Nuclei));
    
    phi0_Unw        = S_parfor.Value{ii}.phi0N(Unwanted_Nuclei);
    phi1_Unw        = S_parfor.Value{ii}.phi1N(Unwanted_Nuclei);
    n0_Unw          = S_parfor.Value{ii}.n0(Unwanted_Nuclei,:);
    n1_Unw          = S_parfor.Value{ii}.n1(Unwanted_Nuclei,:);

    temp            = SubClass_Ent_and_Fid;
    Infid(ii)       = temp.Gate_Infid(length(Unwanted_Nuclei)+GHZsize-1,...
                                      GHZsize-1,phi0_Unw,phi1_Unw,n0_Unw,n1_Unw);
                 
end

clear S_parfor;

to_remove=[];
for Cases = 1:Final_Cases
    
    S{Cases}.Infid = Infid(Cases);
    
    if Infid(Cases).Infid>Infid_Tol
       
        to_remove=[to_remove,Cases];
        
    end
    
end

S(to_remove)=[];

if length(S)<Final_Cases
   
    disp('Excluded additional cases due to gate error tolerances.')
    
end
disp('====================================================================')
disp(['Final GHZ instances:',num2str(length(S))])
disp('====================================================================')

if isempty(S)
    
    warning(['Did not find acceptable case for GHZsize=',num2str(GHZsize),' using the multispin protocol.'])
    
end


end

