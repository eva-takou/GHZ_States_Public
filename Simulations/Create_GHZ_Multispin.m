function S=Create_GHZ_Multispin(GHZsize)
%Script to create GHZM states with the multispin protocol.
%GHZsize: Size of the GHZ state (>=4).

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

cnt=0;
for ii=1:length(OUT)
    
    if isstruct(OUT{ii}) %Not nan-found some acceptable cases.
        cnt=cnt+1;
        TEMP{cnt}=OUT{ii};
    end
    
end
clear OUT

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
for comb=1:L
    
    cnt=0;
    Spins_To_Test = NuclearCombs(comb,:);
    
    for Case = 1:Max_Case
    
        if OUT_new.Target_Nuclei(Case,:)==Spins_To_Test
            
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

        M_Tangle(jj) = prod(S{comb}.Ep(jj,S{comb}.Target_Nuclei));

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
        n0n1(Spin_Indx) = dot(NucSpin.axes{1},NucSpin.axes{2});
        
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

function OUT=Get_ep_for_fixed_k(resonance,wL,s0,s1,A,B,Time_Max,...
                       Unwanted_One_Tangle_Tol,Target_One_Tangle_Tol,GHZsize,Length_Of_Time_Vals)

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

disp(['times I will loop for:',num2str(Lt)])

tic
parfor ii=1:Lt
    
    for Spin_Indx=1:Nnuc
        
       numSpin            = 1;    
       NucSpin            = SubClass_U4Operations(wL,A(Spin_Indx),B(Spin_Indx),s0,s1,numSpin,resonance,1);    
       NucSpin            = NucSpin.CPMG(tt(ii));      
       NucSpin            = NucSpin.Rot_Angles_And_Axes;
       phi(ii,Spin_Indx)  = NucSpin.angles{1};
       n0n1(ii,Spin_Indx) = dot(NucSpin.axes{1},NucSpin.axes{2});
        
    end
    
end

%This has been checked and is correct. Agrees with Makhlin_Inv approach
parfor ii=1:Lt
    
    for NN=1:NMAX
     
        phiN        = acos(cos(NN.*phi(ii,:)));
        G1          = ( cos(phiN/2).^2 + n0n1(ii,:).*sin(phiN/2).^2 ).^2;
        Ep(NN,ii,:) = 1-G1; 
     
    end
end
toc      
clear phi; clear n0n1;
Ep = reshape(Ep,[length(Iters),Nnuc]);

%for some N we exceeded the gate time
indices =  Times.*Iters<Time_Max;
Times   =  Times(indices);
Iters   =  Iters(indices);
Ep      =  Ep(indices,:);

Max_Case  = length(Times);
cnt       = 0;

disp('Checking acceptable cases based on one-tangles tolerances...')
parfor Case = 1 : Max_Case
   
    ep_All    = Ep(Case,:);
    ep_T      = ep_All(ep_All>Target_One_Tangle_Tol);
    ep_U      = ep_All(ep_All<Target_One_Tangle_Tol);
 
    if length(ep_T)==GHZsize-1 && isempty(ep_U(ep_U>Unwanted_One_Tangle_Tol)) 
                               
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

