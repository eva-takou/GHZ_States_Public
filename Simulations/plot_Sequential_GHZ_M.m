function plot_Sequential_GHZ_M
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
FNTsize=16;
legs={'GHZ_4','GHZ_5','GHZ_6','GHZ_7','GHZ_8','GHZ_9','GHZ_{10}'};

%========================================================================
close all;


load('/Users/evatakou/Documents/MATLAB/Defect_Nuclear_GHZ_states/GHZ_states_Preparation/Simulations/GHZ_Data_of_Sequential_NEW/GHZ4_Sequential.mat');
S{1}=OUT;

load('/Users/evatakou/Documents/MATLAB/Defect_Nuclear_GHZ_states/GHZ_states_Preparation/Simulations/GHZ_Data_of_Sequential_NEW/GHZ5_Sequential.mat');
S{2}=OUT;

load('/Users/evatakou/Documents/MATLAB/Defect_Nuclear_GHZ_states/GHZ_states_Preparation/Simulations/GHZ_Data_of_Sequential_NEW/GHZ6_Sequential.mat');
S{3}=OUT;

load('/Users/evatakou/Documents/MATLAB/Defect_Nuclear_GHZ_states/GHZ_states_Preparation/Simulations/GHZ_Data_of_Sequential_NEW/GHZ7_Sequential.mat');
S{4}=OUT;

load('/Users/evatakou/Documents/MATLAB/Defect_Nuclear_GHZ_states/GHZ_states_Preparation/Simulations/GHZ_Data_of_Sequential_NEW/GHZ8_Sequential.mat');
S{5}=OUT;

load('/Users/evatakou/Documents/MATLAB/Defect_Nuclear_GHZ_states/GHZ_states_Preparation/Simulations/GHZ_Data_of_Sequential_NEW/GHZ9_Sequential.mat');
S{6}=OUT;

load('/Users/evatakou/Documents/MATLAB/Defect_Nuclear_GHZ_states/GHZ_states_Preparation/Simulations/GHZ_Data_of_Sequential_NEW/GHZ10_Sequential.mat');
S{7}=OUT;


for ii=1:length(S)
    
[One_Tangles{ii},M_Way_Ep{ii},Gate_Time{ii},Gate_Error{ii}]=Get_Data(S{ii});    
    
end

%=========================================================================
rows = 7;
cols = 1;
figure(1)

%======== Plot of one-tangles ===========================================
for ii=1:rows
     
    subplot(rows,cols,ii) 
    
    b=bar([1:length(S{ii}.Target_Nuclei)],One_Tangles{ii},...
        'BarWidth',0.5,'linewidth',.8);
    
    for jj=1:length(b)
        b(jj).FaceColor = Magentas{1};
    end
    
    ax=gca;
    
    max_case  = length(S{ii}.Target_Nuclei);
    
    for kk=1:max_case
        
        Lnuc = length(S{ii}.Target_Nuclei{kk});
        
            for jj=1:Lnuc
                
                
                if ii==1 %GHZ4
                    
                    fntSizeCarbon = 11;
                    
                elseif ii==2 %GHZ5
                    
                    fntSizeCarbon = 11;
                    
                elseif ii==3 %GHZ6
                    
                    fntSizeCarbon = 11;
                    
                elseif ii==4 %GHZ7
                    
                    fntSizeCarbon = 8;
                    
                elseif ii==5  %GHZ8
                    
                    fntSizeCarbon = 8;
                    
                elseif ii==6 %GHZ9
                    
                    fntSizeCarbon = 10;
                   
                elseif ii==7 %GHZ10
                    
                    fntSizeCarbon = 12;
                    
                end
                    
                
             if (ii==4 && (-1)^jj==1) || (ii==5 && (-1)^jj==1)
                 
                text(ax,'String',strcat("C",num2str(S{ii}.Target_Nuclei{kk}(jj))),...
                 'Position',[b(jj).XEndPoints(kk),max(ax.YLim)+0.025],...
                 'color','k','fontsize',fntSizeCarbon,...
                 'Rotation',90,'fontname','microsoft sans serif'); %'interpret','latex'
             
                 
             elseif (ii==4 && (-1)^jj==-1) || (ii==5 && (-1)^jj==-1)
                 
                text(ax,'String',strcat("C",num2str(S{ii}.Target_Nuclei{kk}(jj))),...
                 'Position',[b(jj).XEndPoints(kk),max(ax.YLim)+0.002],...
                 'color','k','fontsize',fntSizeCarbon,...
                 'Rotation',90,'fontname','microsoft sans serif'); %'interpret','latex'
                 
             else 
                 
                text(ax,'String',strcat("C",num2str(S{ii}.Target_Nuclei{kk}(jj))),...
                 'Position',[b(jj).XEndPoints(kk),max(ax.YLim)+0.002],...
                 'color','k','fontsize',fntSizeCarbon,...
                 'Rotation',90,'fontname','microsoft sans serif'); %'interpret','latex'
                 
             end
             
             
            end
            
        
    end
    
    ylim([min(One_Tangles{ii},[],'all')-0.02,1])
    
    fig_props(FNTsize,'','$\epsilon_p^{nuc.}$')
    title(legs{ii},'fontsize',15)
end
xlabel('Case #')

rows = 7;
cols = 3;

%================ Fig 2 ===================================================
figure(2)

%======== Plot of M-way ep ===============================================
for ii=1:rows
    
    subplot(rows,cols,1+(ii-1)*cols)  %indx: 2, 5, 8, 12,...
    
    b=bar([1:length(S{ii}.Target_Nuclei)],M_Way_Ep{ii},...
        'BarWidth',0.5,'linewidth',.8);
    
    for jj=1:length(b)
        b(jj).FaceColor = Blues{1};
    end
    
    ylim([min(M_Way_Ep{ii},[],'all')-0.02,1])
    fig_props(FNTsize,'Case #','$\epsilon_{p,M}(U)$')
    
    
    title(legs{ii},'fontsize',15)
    
end


%======== Plot of gate errors ===========================================

for ii=1:rows
    
    subplot(rows,cols,2+(ii-1)*cols)
    
    b=bar([1:length(S{ii}.Target_Nuclei)],Gate_Error{ii},...
        'BarWidth',0.5,'linewidth',.8);
    
    for jj=1:length(b)
        b(jj).FaceColor = gate_err_col;
    end
    
    fig_props(FNTsize,'Case #','$1-F$')
    ylim([0,max(Gate_Error{ii},[],'all')+0.005])
    
end

%======== Plot of gate times ===========================================

for ii=1:rows
    
    subplot(rows,cols,3+(ii-1)*cols)
    
    b=bar([1:length(S{ii}.Target_Nuclei)],Gate_Time{ii}/1000,...
        'BarWidth',0.5,'linewidth',.8);
    
    for jj=1:length(b)
        b(jj).FaceColor = gate_time_col;
    end
    
    fig_props(FNTsize,'Case #','$T$ (ms)')
    
end


end

function [One_Tangles,M_Way_Ep,Gate_Time,Gate_Error]=Get_Data(input_struct)

S=input_struct;


for ii=1:length(S.Target_Nuclei)
   
    for kk=1:length(S.EP_Target{ii})
        
       One_Tangles(ii,kk)= S.EP_Target{ii}(kk);
        
    end
    
    Rot_Angles(ii,:) = S.phi0N{ii}(S.Target_Nuclei{ii});
    N0N1(ii,:)       = S.n0n1{ii}(S.Target_Nuclei{ii});
    M_Way_Ep(ii)     = prod(S.EP_Target{ii});
    Gate_Time(ii)    = S.Total_Time{ii};
    Gate_Error(ii)   = S.Infid{ii}.Infid;
    
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