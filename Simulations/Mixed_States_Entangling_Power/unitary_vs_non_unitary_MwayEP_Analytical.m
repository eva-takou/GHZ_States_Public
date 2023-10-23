function unitary_vs_non_unitary_MwayEP_Analytical(path_Sequential,path_Multispin)
%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
%
%Script to compare the unitary with the non-unitary M-way entangling power.
%Input: path of sequential or multispin data.
%       The data should have all the info about target nuclei, unwanted
%       nuclei, iterations/times of sequence.
%
%       This script uses instead of numerical the analytical closed-form
%       expression for the non-unitary M-tangling power.

[ep_U_Seq,ep_E_Seq,ep_E_Aprox_Seq]=unpack_data('Sequential',path_Sequential);
[ep_U_Mul,ep_E_Mul,ep_E_Aprox_Mul]=unpack_data('Multispin',path_Multispin);

%=============== Plot the results =========================================
close all;  clc;

Magentas = {[0.4940 0.1840 0.5560],[147,112,219]/255};
Blues    = {[0,104,139]/255,[25,25,112]/255,[65,105,225]/255,[39,64,139]/255,[0,0,128]/255};
Green    = {[60,179,113]/255};
rows    = 2;
cols    = 1;
fntsize = 20;

%========= Plot the sequential first =====================================

subplot(rows,cols,1)

b=bar(1:length(ep_U_Seq),[ep_U_Seq;ep_E_Seq;ep_E_Aprox_Seq]);
b(1).FaceColor=Blues{1};
%b(1).FaceAlpha=0.2;
b(2).FaceColor=Magentas{1};
b(2).FaceAlpha=0.4;
b(3).FaceColor=Green{1};
b(3).FaceAlpha=0.4;


% hold on
% b=bar(1:length(ep_E_Seq),ep_E_Seq);
% b.FaceColor=Blues{1};
ylim([0.8,1])
xlab='Case #';
ylab='$\epsilon_{p,M}$';
legend({'$\epsilon_{p,M}(U)$','$\epsilon_{p,M}(\mathcal{E})$','$\epsilon_{p,M}(\mathcal{E})$ Approx.'},...
    'interpreter','latex','location','best','color','none',...
    'edgecolor','none','NumColumns',3)


% hold on 
% b=bar(1:length(ep_E_Aprox_Seq),ep_E_Aprox_Seq);
% b.FaceColor='g';
% b.FaceAlpha=0.1;
% 
fig_props(fntsize,xlab,ylab)


%========== Plot the multispin =========================================


subplot(rows,cols,2)

b=bar(1:length(ep_U_Mul),[ep_U_Mul;ep_E_Mul;ep_E_Aprox_Mul]);
b(1).FaceColor=Blues{1};
b(2).FaceColor=Magentas{1};
b(2).FaceAlpha=0.4;
b(3).FaceColor=Green{1};
b(3).FaceAlpha=0.4;
hold on
%b=bar(1:length(ep_E_Mul),ep_E_Mul);
%b.FaceColor=Blues{1};
%ylim([0.75,1])

xlab='Case #';
ylab='$\epsilon_{p,M}$';


% hold on 
% b=bar(1:length(ep_E_Aprox_Mul),ep_E_Aprox_Mul);
% b.FaceColor='g';
% b.FaceAlpha=0.1;
fig_props(fntsize,xlab,ylab)


end



function [ep_U,ep_E,ep_E_Aprox]=unpack_data(Scheme,path)

M = 4;  
d = 2;  
load(path,"OUT");
S = OUT;

tol = 4e-3;

switch Scheme
    
    case 'Sequential'
        
        Max_Case = length(S.Target_Nuclei);
        
    parfor ii=1:Max_Case

        disp(['running case=',num2str(ii),' out of ',num2str(Max_Case)])
        
        Target_Nuclei   = S.Target_Nuclei{ii};
        Unwanted_Nuclei = S.Unwanted_Nuclei_Names{ii};
        temp            = S.EP_Unwanted{ii};
        indices         = temp>tol; 
        Unwanted_Nuclei = Unwanted_Nuclei(indices);
        Ntarget         = length(Target_Nuclei);
        Nnuc            = length(Target_Nuclei)+length(Unwanted_Nuclei);

        %Parameters of unwanted spins evolution
        phi0 = S.phi0N{ii}(Unwanted_Nuclei);
        phi1 = S.phi1N{ii}(Unwanted_Nuclei);
        n0   = S.n0{ii}(Unwanted_Nuclei,:);
        n1   = S.n1{ii}(Unwanted_Nuclei,:);

        %Parameters of target spins evolution
        th0  = S.phi0N{ii}(Target_Nuclei);
        th1  = S.phi1N{ii}(Target_Nuclei);
        u0   = S.n0{ii}(Target_Nuclei,:);
        u1   = S.n1{ii}(Target_Nuclei,:);

        temp     = SubClass_Ent_and_Fid;
        temp     = temp.NonUniMwayEp_CR_Analytical(Nnuc,Ntarget,phi0,phi1,n0,n1,th0,th1,u0,u1);
        ep_E(ii) = temp.epM_nonuni/(d/(d+1))^M; %Scale by max value
        ep_U(ii) = prod(S.EP_Target{ii}); 

        temp           = temp.NonUniMwayEp_CR_Aproximate(Nnuc,Ntarget,phi0,phi1,n0,n1,th0,th1,u0,u1);
        ep_E_Aprox(ii) = temp.epM_nonuni/(d/(d+1))^M;
    end
        

        
    case 'Multispin'
        
        Max_Case = length(S);
        
        parfor ii=1:Max_Case

            Target_Nuclei   = S{ii}.Target_Nuclei;
            Unwanted_Nuclei = S{ii}.Unwanted_Nuclei;
            temp            = S{ii}.Ep_Unwanted;
            indices         = temp>tol;
            Unwanted_Nuclei = Unwanted_Nuclei(indices);

            Ntarget         = length(Target_Nuclei);
            Nnuc            = length(Target_Nuclei)+length(Unwanted_Nuclei);

            phi0 = S{ii}.phi0N(Unwanted_Nuclei);
            phi1 = S{ii}.phi1N(Unwanted_Nuclei);

            n0   = S{ii}.n0(Unwanted_Nuclei,:);
            n1   = S{ii}.n1(Unwanted_Nuclei,:);

            th0  = S{ii}.phi0N(Target_Nuclei);
            th1  = S{ii}.phi1N(Target_Nuclei);

            u0   = S{ii}.n0(Target_Nuclei,:);
            u1   = S{ii}.n1(Target_Nuclei,:);

            temp     = SubClass_Ent_and_Fid;
            temp     = temp.NonUniMwayEp_CR_Analytical(Nnuc,Ntarget,phi0,phi1,n0,n1,th0,th1,u0,u1);
            ep_E(ii) = temp.epM_nonuni/(d/(d+1))^M;
            ep_U(ii) = prod(S{ii}.Ep_Target);

            temp           = temp.NonUniMwayEp_CR_Aproximate(Nnuc,Ntarget,phi0,phi1,n0,n1,th0,th1,u0,u1);
            ep_E_Aprox(ii) = temp.epM_nonuni/(d/(d+1))^M;
            
            
        end   

        
end




end



