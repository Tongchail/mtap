%% Plot the cmp points to Opx-Cpx-An rock diagram to define the cmp name.

addpath('../src')

run(['../cal/cal_',calID]);  % load melt model calibration

cal.cmp_mem_pro = cal.cmp_mem(1:end-2,:);

ant = cal.cmp_mem_pro(:,1);
alb = cal.cmp_mem_pro(:,2);
san = cal.cmp_mem_pro(:,3);

chy = cal.cmp_mem_pro(:,7);
fhy = cal.cmp_mem_pro(:,8);
hyp = cal.cmp_mem_pro(:,9);

mau = cal.cmp_mem_pro(:,10);
fau = cal.cmp_mem_pro(:,11);
aug = cal.cmp_mem_pro(:,12);

An  = (ant + alb + san) ./ sum(ant + alb + san + mau + fau + aug + chy + fhy + hyp);
Cpx = (mau + fau + aug) ./ sum(ant + alb + san + mau + fau + aug + chy + fhy + hyp);
Opx = (chy + fhy + hyp) ./ sum(ant + alb + san + mau + fau + aug + chy + fhy + hyp);

figure(30); clf
OpxCpxAn; hold on

[xp, yp] = tern(Cpx, An, Opx); 

idx = 1:3;   % set which points to plot

scatter(xp(idx), yp(idx), 80, 'filled', 'MarkerEdgeColor', 'k'); 

labels = {'cmp1','cmp2','cmp3','cmp4','cmp5'};  

for i = idx
    text(xp(i)+0.02, yp(i)-0.01, labels{i}, ...
         'FontSize', 10, 'Color', 'b', ...
         'HorizontalAlignment','left', ...
         'VerticalAlignment','bottom');
end