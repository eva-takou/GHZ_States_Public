function S=Create_GHZ3_Multispin

%============= Set the basic parameters ====================================
[s0,s1,wL] = load_fixed_params;  % Electron's spin projections, and nuclear larmor freq
[A,B]      = Get_HF_Delft;       % HF parameters of 27 spins
Nnuc       = length(A);          % Number of nuclei
GHZsize    = 3;

[Target_One_Tangle_Tol,Unwanted_One_Tangle_Tol,Infid_Tol,Time_Max,~,Time_Vals]=Tolerances_Multispin(GHZsize);

%=============== Run for multiple resonances ==============================

OUT{1}  = Get_ep_for_fixed_k(1,wL,s0,s1,A,B,Time_Max,Time_Vals);
disp('Finished with 1st resonance.')
OUT{2}  = Get_ep_for_fixed_k(2,wL,s0,s1,A,B,Time_Max,Time_Vals);
disp('Finished with 2nd resonance.')
OUT{3}  = Get_ep_for_fixed_k(3,wL,s0,s1,A,B,Time_Max,Time_Vals);
disp('Finished with 3rd resonance.')
OUT{4}  = Get_ep_for_fixed_k(4,wL,s0,s1,A,B,Time_Max,Time_Vals);
disp('Finished with 4th resonance.')
OUT{5}  = Get_ep_for_fixed_k(5,wL,s0,s1,A,B,Time_Max,Time_Vals);
disp('Finished with 5th resonance.')
OUT{6}  = Get_ep_for_fixed_k(6,wL,s0,s1,A,B,Time_Max,Time_Vals);
disp('Finished with 6th resonance.')
OUT{7}  = Get_ep_for_fixed_k(7,wL,s0,s1,A,B,Time_Max,Time_Vals);
disp('Finished with 7th resonance.')
OUT{8}  = Get_ep_for_fixed_k(8,wL,s0,s1,A,B,Time_Max,Time_Vals);
disp('Finished with 8th resonance.')
OUT{9}  = Get_ep_for_fixed_k(9,wL,s0,s1,A,B,Time_Max,Time_Vals);
disp('Finished with 9th resonance.')
OUT{10} = Get_ep_for_fixed_k(10,wL,s0,s1,A,B,Time_Max,Time_Vals);
disp('Finished with 10th resonance.')

%unfold all cases for each resonance in the following way:
cnt=0;
for Res = 1:length(OUT)
    
    if isstruct(OUT{Res})
        
        Max_Case = length(OUT{Res}.Times);

        for Case = 1:Max_Case
            
            cnt=cnt+1;
            OUT_new.Times(cnt) = OUT{Res}.Times(Case);
            OUT_new.Iters(cnt) = OUT{Res}.Iters(Case);
            OUT_new.Ep(cnt,:)  = OUT{Res}.Ep(Case,:); %size length(cases) x length(spins)

        end
    
    end
end

clear OUT; %to free memory

NuclearCombs = nchoosek(1:24,GHZsize);
[L,~]        = size(NuclearCombs);
Max_Case     = length(OUT_new.Times);

S = parallel.pool.Constant(OUT_new);

%Get the acceptable cases based on tolerances.

parfor Case = 1:Max_Case
    
    ep_All = S.Value.Ep(Case,:);
    ep_T   = ep_All(ep_All>Target_One_Tangle_Tol);
    ep_U   = ep_All(ep_All<Unwanted_One_Tangle_Tol);
    
    if length(ep_T)==GHZsize && length(ep_U)==27-GHZsize
        
        OUT_new2{Case}.Target_Nuclei = find(ep_All>Target_One_Tangle_Tol);
        OUT_new2{Case}.Ep            = S.Value.Ep(Case,:)
        OUT_new2{Case}.Times         = S.Value.Times(Case);
        OUT_new2{Case}.Iters         = S.Value.Iters(Case);
    end
    
end

S = cell(1,L);

%Arange the acceptable cases based on unique nuclear spin combinations.
disp('Aranging acceptable cases based on unique nuclear combinations...')
for comb=1:L

    cnt          = 0;
    Spins_To_Test = NuclearCombs(comb,:);
    
    for Case = 1:length(OUT_new2)
    
        if isstruct(OUT_new2{Case}) & OUT_new2{Case}.Target_Nuclei==Spins_To_Test
            
            cnt=cnt+1;
            S{comb}.Target_Nuclei      = Spins_To_Test;
            S{comb}.Times(cnt)         = OUT_new2{Case}.Times;
            S{comb}.Iters(cnt)         = OUT_new2{Case}.Iters;
            S{comb}.Ep(cnt,:)          = OUT_new2{Case}.Ep;
        end
    end
    
end
disp('Done.')
clear OUT_new2;


to_remove=[];
for ii=1:length(S)
    
    if ~isstruct(S{ii})
       
        to_remove=[to_remove,ii];
        
    end
    
end

S(to_remove)=[];

disp('===================================================================')
disp(['Total unique spin combinations:',num2str(length(S))])
disp('===================================================================')

%For each unique combination select the optimal case (maximum M-tangle)

for comb=1:length(S)
        
    Max_Case  = length(S{comb}.Times);
    M_Tangle  = zeros(1,Max_Case);
    
    for jj=1:Max_Case

        M_Tangle(jj) = prod(S{comb}.Ep(jj,S{comb}.Target_Nuclei));

    end
    
    [~,indx] = max(M_Tangle,[],'all','linear');

    opt_t = S{comb}.Times(indx);
    opt_N = S{comb}.Iters(indx);
    
    S{comb}.Unwanted_Nuclei = setxor(1:27,S{comb}.Target_Nuclei);
    S{comb}.Times           = opt_t;
    S{comb}.Iters           = opt_N;
    S{comb}.Total_Time      = opt_t*opt_N;

    S{comb}.Ep          = S{comb}.Ep(indx,:);
    S{comb}.Ep_Target   = S{comb}.Ep(S{comb}.Target_Nuclei);
    S{comb}.Ep_Unwanted = S{comb}.Ep(S{comb}.Unwanted_Nuclei);

    S{comb}.A_Target    = A(S{comb}.Target_Nuclei);
    S{comb}.B_Target    = B(S{comb}.Target_Nuclei);
    S{comb}.A_Unwanted  = A(S{comb}.Unwanted_Nuclei);
    S{comb}.B_Unwanted  = B(S{comb}.Unwanted_Nuclei);
          
end

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
        NucSpin   = SubClass_U4Operations(wL,A(Spin_Indx),B(Spin_Indx),s0,s1,numSpin,1,S{ii}.Iters);                               
        NucSpin   = NucSpin.CPMG(S{ii}.Times);
        NucSpin   = NucSpin.Rot_Angles_And_Axes;

        phi0N(Spin_Indx) = NucSpin.angles{1};
        phi1N(Spin_Indx) = NucSpin.angles{2};
        n0(Spin_Indx,:)  = NucSpin.axes{1};
        n1(Spin_Indx,:)  = NucSpin.axes{2};
        n0n1(Spin_Indx)  = dot(NucSpin.axes{1},NucSpin.axes{2});
        
    end
    
        S{ii}.phi0N = phi0N;
        S{ii}.phi1N = phi1N;
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
    
    phi0_Unw = S_parfor.Value{ii}.phi0N(Unwanted_Nuclei);
    phi1_Unw = S_parfor.Value{ii}.phi1N(Unwanted_Nuclei);
    n0_Unw   = S_parfor.Value{ii}.n0(Unwanted_Nuclei,:);
    n1_Unw   = S_parfor.Value{ii}.n1(Unwanted_Nuclei,:);

    temp      = SubClass_Ent_and_Fid;
    Infid(ii) = temp.Analy_Expr_Fid2(length(Unwanted_Nuclei)+GHZsize-1,...
                     GHZsize-1,phi0_Unw,phi1_Unw,n0_Unw,n1_Unw);
                 
end
clear S_parfor;

to_remove=[];

for Cases=1:Final_Cases
    
    S{Cases}.Infid=Infid(Cases);
    
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


end


function OUT = Get_ep_for_fixed_k(resonance,wL,s0,s1,A,B,Time_Max,Time_Vals)

Nnuc       = length(A);
unit_times = zeros(1,Nnuc);

for fixed_Nuc = 1:24

    unit_times(fixed_Nuc) = Resonance_Time(resonance,fixed_Nuc,wL,A,B,s0,s1);
    
end

tmin                 = (1-0.25)    * min(unit_times); 
tmax                 = (1+0.25)    * max(unit_times); 
tt                   = linspace(tmin,tmax,Time_Vals); 
Length_Of_Time_Vals  = length(tt);
Number_Of_Maxima     = 50; %How many maxima of one-tangles to keep later.

%For each time, get the maxima of each nuclear spin. Then, 
%collect all the iterations as long as N*tt<Time_Max
%Then, use all those N to get ep.

cnt=0;
for kk=1:Length_Of_Time_Vals
    
    for Spin_Indx = 1:24
        
        numSpin = 1;    
        NucSpin = SubClass_U4Operations(wL,A(Spin_Indx),B(Spin_Indx),s0,s1,numSpin,resonance,1);                               
        NucSpin = NucSpin.CPMG(tt(kk));
        NucSpin = NucSpin.Rot_Angles_And_Axes;
        
        if dot(NucSpin.axes{1},NucSpin.axes{2})>=0
            continue
        end
        
        NucSpin = NucSpin.Expected_Maxima(Number_Of_Maxima);
        temp_N  = NucSpin.Nmax;
        
        for ii=1:length(temp_N)
            
            if temp_N(ii)*tt(kk)<Time_Max
                
                cnt        = cnt+1;
                Times(cnt) = tt(kk);                 
                Iters(cnt) = temp_N(ii);
                
            end
            
        end
        
    end

end


if cnt==0 %we failed
    
    OUT=nan;
    
    return
end

Max_Cases = cnt;

AHF = parallel.pool.Constant(A);
BHF = parallel.pool.Constant(B);

parfor Case = 1 : Max_Cases
    
    for Spin_Indx = 1:Nnuc
        
      numSpin = 1;    
      NucSpin = SubClass_U4Operations(wL,AHF.Value(Spin_Indx),BHF.Value(Spin_Indx),...
                                      s0,s1,numSpin,resonance,Iters(Case));  
      NucSpin                 = NucSpin.CPMG(Times(Case));
      NucSpin                 = NucSpin.Makhlin_Inv;
      Ep(Case,Spin_Indx)      = NucSpin.Ep/(2/9);
        
    end
    
end


OUT.Times=Times;
OUT.Iters=Iters;
OUT.Ep = Ep;


end



