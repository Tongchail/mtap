%% Phase fraction of [Mineral] vs Temperature


%% Plot all phases in one figure but subplots

nphs = size(PHS_frc,2);

figure (1)
for i = 1:nphs
    subplot(ceil(nphs/2),2,i);
    scatter(Tmp, PHS_frc(:,i), 36, 'filled');      % point size 36
    set(gca, 'XDir','reverse')  
    xlabel('Temperature (°C)');
    ylabel('Phase fraction (wt%)');
    title(char(phs(i)));
    grid on;
end
sgtitle('Phase fractions vs Temperature');


%%  different colour
load MtAp_calibration_all.mat
nphs = size(PHS_frc,2);
N    = size(PHS_frc,1);

pattern = [34, 32, 31, 34, 34];   
K = numel(pattern);
cmap = lines(K);             

figure(1)
for i = 1:nphs
    subplot(ceil(nphs/2),2,i); hold on

    idx_start = 1;
    shown = false(1,K);        
    seg = 1;                    

    while idx_start <= N
        k   = mod(seg-1, K) + 1;               
        len = pattern(k);
        idx_end = min(idx_start + len - 1, N);  
        if ~shown(k)           
            scatter(Tmp(idx_start:idx_end), PHS_frc(idx_start:idx_end,i), ...
                    36, cmap(k,:), 'filled', 'DisplayName', sprintf('Group %d', k));
            shown(k) = true;
        else        
            scatter(Tmp(idx_start:idx_end), PHS_frc(idx_start:idx_end,i), ...
                    36, cmap(k,:), 'filled', 'HandleVisibility', 'off');
        end
        idx_start = idx_end + 1;
        seg = seg + 1;
    end

    set(gca, 'XDir','reverse')
    xlabel('Temperature (°C)');
    ylabel('Phase fraction (wt%)');

    if exist('phs','var')
        if iscell(phs), title(phs{i});
        else,           title(string(phs(i)));
        end
    else
        title(sprintf('Phase %d', i));
    end

    grid on
    legend('show','Location','best'); 
    hold off
end

sgtitle('Phase fractions vs Temperature');  


%% remove liq  % plot all in one figure

mask = ~(phs == "liq"); 
PHS_frc_new = PHS_frc(:,mask);
phs_new = phs(mask);

figure (2);
hold on;

nphs = size(PHS_frc_new,2);
for i = 1:nphs
    scatter(Tmp, PHS_frc_new(:,i), 36, 'filled', 'DisplayName', char(phs_new(i)));
end

set(gca, 'XDir','reverse');  
xlabel('Temperature (°C)');
ylabel('Phase fraction (%)');
title('Phase fractions vs Temperature (excluding liq and qtz)');
legend('show','Location','best');
grid on;
hold off;

%% remove liq and q   in one figure but subplots

mask = ~(phs == "liq"); 
PHS_frc_new = PHS_frc(:,mask);
phs_new = phs(mask);

nphs = size(PHS_frc_new,2);
nrow = ceil(sqrt(nphs));   % row
ncol = ceil(nphs/nrow);    % column

figure (3);
for i = 1:nphs
    subplot(nrow, ncol, i);
    scatter(Tmp, PHS_frc_new(:,i), 36, 'filled');
    set(gca, 'XDir','reverse');   
    xlabel('T (°C)');
    ylabel('Fraction (%)');
    title(char(phs_new(i)));
    grid on;
end
sgtitle('Phase fractions vs Temperature (excluding liq & qtz)');