function plot_MultiSpin_GHZ3(S)
%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
%
%Function to plot the optimal GHZ3 cases obtained based on the methods developed in
%arXiv:2302.05580, for the multispin scheme.
%--------------------------------------------------------------------------
%Input: S: a struct with various fields that stores the optimal cases and
%          information about the total evolution.
%--------------------------------------------------------------------------

close all;

%========== Colors ========================================================

gate_err_col  = [0.4,0,0.2];
gate_time_col = [0.4,0.7,1];

Blues={[0,104,139]/255,[25,25,112]/255,[65,105,225]/255,[39,64,139]/255,[0,0,128]/255};

Magentas={[0.4940 0.1840 0.5560],[147,112,219]/255};

Yellows={[255,185,15]/255,[255,215,0]/255,[238,238,0]/255};

Purples={[204,153,255]/255,[153,0,153]/255,[102,0,51]/255,...
        [51,0,102]/255};

Greens={[127,255,212]/255,[102,205,170]/255,[69,139,116]/255}  ;

Grays={[119,136,153]/255,[159,182,205]/255,[112,128,144]/255};

Reds={[205,55,0]/255,[205,38,38]/255,[139,0,0]/255};

OtherColors={[205,179,139]/255};

clrs_Errors={[204,255,255]./255,...
             [102,178,255]./255,...
             [0,128,255]./255,...
             [0,76,153]./255,...
             [102,255,102]./255,...
             [0,153,0]./255,...
             [255,153,153]./255,...
             [204,0,0]./255,...
             [102,0,0]./255};


%========= Default options ===============================================


FNTsize=19;
rows = 6;
cols = 1;
%========================================================================


[One_Tangles,Rot_Angles,N0N1,M_Way_Ep,Gate_Time,Gate_Error]=Get_Data(S);


figure(1)

%========== Plot of one-tangles ===========================================
subplot(rows,cols,1)

b=bar([1:length(S)],One_Tangles,'BarWidth',0.5,'linewidth',.8);
for ii=1:length(b)
    
   b(ii).FaceColor = Magentas{1};
   hold on
    
end

ax=gca;
max_case  = length(S);

for kk=1:max_case
    
        Lnuc = length(S{kk}.Target_Nuclei);
        
            for jj=1:Lnuc
                
                text(ax,'String',strcat("C",num2str(S{kk}.Target_Nuclei(jj))),...
                 'Position',[b(jj).XEndPoints(kk),max(ax.YLim)+0.002],...
                 'color','k','fontsize',14,...
                 'Rotation',90,'fontname','microsoft sans serif'); %'interpret','latex'
                
            end

end

ylim([0.9,1]);    xlim([0,length(S)+1]);
fig_props(FNTsize,'Case #','$\epsilon_p^{nuc.}$')


%============= Plot of rot angles =======================================

subplot(rows,cols,2)
b=bar([1:length(S)],Rot_Angles/(pi/2),'BarWidth',0.5,'linewidth',.8);

for ii=1:2
    b(ii).FaceColor=Greens{2};    
end

xlim([0,length(S)+1]); %ylim([0.8,1])

fig_props(FNTsize,'Case #','$\phi/(\pi/2)$')

%============= Plot of dot Prod ==========================================

subplot(rows,cols,3)
b=bar([1:length(S)],N0N1,'BarWidth',0.5,'linewidth',.8);
xlim([0,length(S)+1]); ylim([-1,0.2])

for ii=1:2
    b(ii).FaceColor=Grays{1};
   
end

fig_props(FNTsize,'Case #','$\textbf{n}_0\cdot \textbf{n}_1$')

%=============== Plot of M-way Ep =========================================

subplot(rows,cols,4)

b=bar([1:length(S)],M_Way_Ep,'BarWidth',0.5,'linewidth',0.8);

b.FaceColor=Blues{1};

fig_props(FNTsize,'Case #','$\epsilon_{p,M}(U)$')
ylim([0.84,1]);    xlim([0,length(S)+1]);
set(gca,'ytick',[0.9,1])
%========== Plot of Gate Time =============================================

subplot(rows,cols,5)
b=bar([1:length(S)],Gate_Time/1e3,'BarWidth',0.5,'linewidth',0.8);
b.FaceColor=gate_time_col;
xlim([0,length(S)+1]);
fig_props(FNTsize,'Case #','$T$ (ms)')

%======== Plot of Gate Error ==============================================


subplot(rows,cols,6)
b=bar([1:length(S)],Gate_Error,'stacked','BarWidth',0.5,'linewidth',0.8);
b(1).FaceColor=gate_err_col;

line([0,length(S)+1],[0.05,0.05],'linewidth',0.5,'linestyle','--','color','k')
fig_props(FNTsize,'Case #','$1-F$')
ylim([0,0.11]); xlim([0,length(S)+1]);

%ax=gca;
%set(ax,'YTick',[0,0.05,0.1],'yticklabels',{'0','0.1','0.2'})


end


function [One_Tangles,Rot_Angles,N0N1,M_Way_Ep,Gate_Time,Gate_Error]=Get_Data(input_struct)

S = input_struct;


for ii=1:length(S)
   
    for kk=1:length(S{ii}.Target_Nuclei)
    
         One_Tangles(ii,kk)=S{ii}.Ep_Target(kk);
    
    end
    
    Rot_Angles(ii,:) = S{ii}.phi0N(S{ii}.Target_Nuclei);
    N0N1(ii,:)       = S{ii}.n0n1(S{ii}.Target_Nuclei);
    M_Way_Ep(ii)     = prod(S{ii}.Ep_Target);
    Gate_Time(ii)    = S{ii}.Total_Time;
    Gate_Error(ii)   = S{ii}.Infid.Infid;
    
end




end



function fig_props(fntsize,xlab,ylab)


set(gca,'fontname','microsoft sans serif','fontsize',fntsize)
set(gcf,'color','w')

xlabel(xlab)
ylabel(ylab,'interpreter','latex')

set(gca, 'Box', 'on', 'TickDir', 'out', 'TickLength', [.003 .003], ...
    'XMinorTick', 'off', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
    'LineWidth', 0.5)

end