%% Plot the cmp points to TAS to define the cmp name.
clc; clear;
addpath('../cal')
addpath('../src')
cal_MtAp_750_new;
% prepare for plotting
TX = {'Interpreter','Latex'}; FS = {'FontSize',12};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',10};
UN = {'Units','Centimeters'};
CL = {'Color',[0.0 0.0 0.0],[0.80 0.15 0.10],[0.10 0.15 0.65],[1.0 0.7 0.0]};
LW = {'LineWidth',2};


TAS; colorbar('off'); hold on;

cal.cmp_oxd_pro = cal.cmp_oxd(:, cal.ioxd);
cal.cmp_oxd_pro = cal.cmp_oxd_pro(1:end-2, :);

cSi =  cal.cmp_oxd_pro(:,1)./sum( cal.cmp_oxd_pro(:,1:end-1),2).*100;
cNK = sum( cal.cmp_oxd_pro(:,[7,8]),2)./sum( cal.cmp_oxd_pro(:,1:end-1),2).*100;
labels = {'cmp1','cmp2','cmp3','cmp4','cmp5'};

scatter( cSi(:), cNK(:),120,'filled','o','MarkerEdgeColor',CL{2},LW{1},1.5);
for i = 1:length(cSi)
    text(cSi(i), cNK(i), labels{i}, ...
         'FontSize', 12, 'Color', 'blue', ...
         'HorizontalAlignment','left', ...
         'VerticalAlignment','bottom');
end

xlabel('SiO$_2$ [wt \%]',TX{:},'FontSize',15); ylabel('Na$_2$O + K$_2$O [wt \%]',TX{:},'FontSize',15);



