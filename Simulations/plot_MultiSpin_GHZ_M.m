function plot_MultiSpin_GHZ_M

close all
%======================= Colors =========================================
Magentas={[0.4940 0.1840 0.5560],[147,112,219]/255};

Blues={[0,104,139]/255,[25,25,112]/255,[65,105,225]/255,[39,64,139]/255,[0,0,128]/255};

Reds={[205,55,0]/255,[205,38,38]/255,[139,0,0]/255};

Greens={[127,255,212]/255,[102,205,170]/255,[69,139,116]/255}  ;

Grays={[119,136,153]/255,[159,182,205]/255,[112,128,144]/255};

OtherColors={[205,179,139]/255};

gate_err_col  = [0.4,0,0.2];
gate_time_col = [0.4,0.7,1];

%====================== Default options =================================
FNTsize=19;
legs={'GHZ_4','GHZ_5','GHZ_6','GHZ_7','GHZ_8','GHZ_9','GHZ_{10}'};

%========================================================================
close all;


load('/Users/evatakou/Documents/MATLAB/Defect_Nuclear_GHZ_states/GHZ_states_Preparation/Simulations/GHZ_Data_of_MultiSpin_NEW/GHZ4_Multispin.mat');
S{1}=OUT;

load('/Users/evatakou/Documents/MATLAB/Defect_Nuclear_GHZ_states/GHZ_states_Preparation/Simulations/GHZ_Data_of_MultiSpin_NEW/GHZ5_Multispin.mat');
S{2}=OUT;

load('/Users/evatakou/Documents/MATLAB/Defect_Nuclear_GHZ_states/GHZ_states_Preparation/Simulations/GHZ_Data_of_MultiSpin_NEW/GHZ6_Multispin.mat');
S{3}=OUT;

load('/Users/evatakou/Documents/MATLAB/Defect_Nuclear_GHZ_states/GHZ_states_Preparation/Simulations/GHZ_Data_of_MultiSpin_NEW/GHZ7_Multispin.mat');
S{4}=OUT;


load('/Users/evatakou/Documents/MATLAB/Defect_Nuclear_GHZ_states/GHZ_states_Preparation/Simulations/GHZ_Data_of_MultiSpin_NEW/GHZ8_Multispin.mat');
S{5}=OUT;

load('/Users/evatakou/Documents/MATLAB/Defect_Nuclear_GHZ_states/GHZ_states_Preparation/Simulations/GHZ_Data_of_MultiSpin_NEW/GHZ9_Multispin.mat');
S{6}=OUT;



for ii=1:length(S)
    
[One_Tangles{ii},M_Way_Ep{ii},Gate_Time{ii},Gate_Error{ii}]=Get_Data(S{ii});    
    
end

%=========================================================================
rows = 6;
cols = 4;
figure(1)

%======== Plot of one-tangles ===========================================
for ii=1:rows
     
    subplot(rows,cols,1+(ii-1)*cols) 
    
    b=bar([1:length(S{ii})],One_Tangles{ii},...
        'BarWidth',0.5,'linewidth',.8);
    
    for jj=1:length(b)
        b(jj).FaceColor = Magentas{1};
    end
    
    ax=gca;
    
    max_case  = length(S{ii});
    
    for kk=1:max_case
        Lnuc = length(S{ii}{kk}.Target_Nuclei);
        
            for jj=1:Lnuc
                
                if ii~=1
                
                    fntSzCarbon=13;
                else
                    fntSzCarbon=10;
                    
                end
                    
                text(ax,'String',strcat("C",num2str(S{ii}{kk}.Target_Nuclei(jj))),...
                 'Position',[b(jj).XEndPoints(kk),max(ax.YLim)+0.002],...
                 'color','k','fontsize',fntSzCarbon,...
                 'Rotation',90,'fontname','microsoft sans serif'); %'interpret','latex'
                
            end
           
    end
    
    
    ylim([min(One_Tangles{ii},[],'all')-0.01,1])
    
    fig_props(FNTsize,'','$\epsilon_p^{nuc.}$')
    
end

xlabel('Case #')



%figure(2)

%======== Plot of M-way ep ===============================================
for ii=1:rows
    
    subplot(rows,cols,2+(ii-1)*cols)  %indx: 2, 5, 8, 12,...
    
    b=bar([1:length(S{ii})],M_Way_Ep{ii},...
        'BarWidth',0.5,'linewidth',.8);
    
    for jj=1:length(b)
        b(jj).FaceColor = Blues{1};
    end
    
    ylim([min(M_Way_Ep{ii},[],'all')-0.03,1])
    fig_props(FNTsize,'','$\epsilon_{p,M}(U)$')
    title(legs{ii},'fontsize',17)
    
end
xlabel('Case #')


%======== Plot of gate errors ===========================================

for ii=1:rows
    
    subplot(rows,cols,3+(ii-1)*cols)
    
    b=bar([1:length(S{ii})],Gate_Error{ii},...
        'BarWidth',0.5,'linewidth',.8);
    
    for jj=1:length(b)
        b(jj).FaceColor = gate_err_col;
    end
    
    %ylim([min(M_Way_Ep{ii},[],'all')-0.02,1])
    fig_props(FNTsize,'','$1-F$')
    ylim([0,max(Gate_Error{ii},[],'all')+0.005])
    
end

xlabel('Case #')
%======== Plot of gate times ===========================================

for ii=1:rows
    
    subplot(rows,cols,4+(ii-1)*cols)
    
    b=bar([1:length(S{ii})],Gate_Time{ii}/1000,...
        'BarWidth',0.5,'linewidth',.8);
    
    for jj=1:length(b)
        b(jj).FaceColor = gate_time_col;
    end
    
    ylim([0,max(Gate_Time{ii}/1000,[],'all')+0.1])
    fig_props(FNTsize,'','$T$ (ms)')
    
end
xlabel('Case #')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





end

function [One_Tangles,M_Way_Ep,Gate_Time,Gate_Error]=Get_Data(input_struct)

S=input_struct;


for ii=1:length(S)
   
    for kk=1:length(S{ii}.Ep_Target)
        
       One_Tangles(ii,kk)= S{ii}.Ep_Target(kk);
        
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


set(gca, 'Box', 'on', 'TickDir', 'out', 'TickLength', [.01 .01], ...
    'XMinorTick', 'off', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
    'LineWidth', 0.5)

end