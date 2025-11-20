%% Mafic rocks ternary base map: Opx–Cpx–An
figure(1); clf; hold on;
axis equal off;

% inverted ternary diagram
tern = @(Cpx,An,Opx) deal( ...
    (0.5*An + Opx) ./ (Cpx + An + Opx), ... 
     (sqrt(3)/2) * (Cpx + Opx) ./ (Cpx + An + Opx)); 

% Points location
[xC,yC] = tern(100,0,0);    % Cpx
[xA,yA] = tern(0,100,0);    % An
[xO,yO] = tern(0,0,100);    % Opx

plot([xA xC xO xA],[yA yC yO yA],'k','LineWidth',1.5);

% Cpx contour line：Cpx = 90, 5
[xc90_l,yc90_l] = tern(90,10,0);
[xc90_r,yc90_r] = tern(90,0,10);
plot([xc90_l xc90_r],[yc90_l yc90_r],'k','LineWidth',1);  % Cpx=90

[xc5_l,yc5_l] = tern(5,90,5);
[xc5_r,yc5_r] = tern(5,5,90);
plot([xc5_l xc5_r],[yc5_l yc5_r],'k','LineWidth',1);  % Cpx=5

% middle line Cpx = Opx（ GABBRO / NORITE）
[xm_b,ym_b] = tern(5,90,5);    % bottom on An=90
[xm_t,ym_t] = tern(45,10,45);  % top on An=10
plot([xm_b xm_t],[ym_b ym_t],'k','LineWidth',1);

% An contour line：An = 90, 10
[xa90_l,ya90_l] = tern(10,90,0);
[xa90_r,ya90_r] = tern(0,90,10);
plot([xa90_l xa90_r],[ya90_l ya90_r],'k','LineWidth',1);  % An=90

[xa10_l,ya10_l] = tern(85,10,5);
[xa10_r,ya10_r] = tern(5,10,85);
plot([xa10_l xa10_r],[ya10_l ya10_r],'k','LineWidth',1);  % An=10

% Opx contour line：Opx = 90, 5
[xo90_l,yo90_l] = tern(10,0,90);
[xo90_r,yo90_r] = tern(0,10,90);
plot([xo90_l xo90_r],[yo90_l yo90_r],'k','LineWidth',1);  % Opx=90

[xco5_l,yco5_l] = tern(90,5,5);
[xco5_r,yco5_r] = tern(5,90,5);
plot([xco5_l xco5_r],[yco5_l yco5_r],'k','LineWidth',1);  % Opx=5


text(xC-0.1,yC,'Cpx','FontSize',12,'FontWeight','bold');
text(xO+0.03,yO,'Opx','FontSize',12,'FontWeight','bold');
text(xA,yA-0.02,'An','FontSize',12,'FontWeight','bold');

title('Mafic rocks (Opx–Cpx–An)','FontSize',12,'FontWeight','bold'); hold on;

% rock type name define
% Anorthosite
[C_ano, A_ano, O_ano] = deal(5, 92, 3);  
[xAno, yAno] = tern(C_ano, A_ano, O_ano);
text(xAno, yAno, 'ANORTHOSITE', ...
     'FontSize', 9, 'FontWeight','bold', ...
     'HorizontalAlignment','center', ...
     'VerticalAlignment','middle');

% OPX-Gabbro 
[C_gab, A_gab, O_gab] = deal(60, 40, 15);
[xGab, yGab] = tern(C_gab, A_gab, O_gab);
text(xGab, yGab, 'OPX-GABBRO', ...
     'FontSize', 9, 'FontWeight','bold', ...
     'HorizontalAlignment','center', ...
     'VerticalAlignment','middle');

% CPX-Norite
[C_nor, A_nor, O_nor] = deal(15, 40, 60);
[xNor, yNor] = tern(C_nor, A_nor, O_nor);
text(xNor, yNor, 'CPX-NORITE', ...
     'FontSize', 9, 'FontWeight','bold', ...
     'HorizontalAlignment','center', ...
     'VerticalAlignment','middle');

% Gabbro–Norite
[C_gn, A_gn, O_gn] = deal(40, 20, 40);
[xGN, yGN] = tern(C_gn, A_gn, O_gn);
text(xGN, yGN, 'GABBRO NORITE', ...
     'FontSize', 9, 'FontWeight','bold', ...
     'HorizontalAlignment','center', ...
     'VerticalAlignment','middle');

% Plagioclase-bearing websterite
[C_pw, A_pw, O_pw] = deal(45, 5, 50);
[xPW, yPW] = tern(C_pw, A_pw, O_pw);
text(xPW, yPW, 'PLAGIOCLASE-bearing WEBSTERITE', ...
     'FontSize', 8, 'FontWeight','bold', ...
     'HorizontalAlignment','center', ...
     'VerticalAlignment','middle');

% Gabbro
[C_gabL, A_gabL, O_gabL] = deal(50, 50, 3);  
[xGabL, yGabL] = tern(C_gabL, A_gabL, O_gabL);

text(xGabL, yGabL, 'GABBRO', ...
     'FontSize', 9, 'FontWeight','bold', ...
     'HorizontalAlignment','center', ...
     'VerticalAlignment','middle', 'Rotation', -60);

% NORITE
[C_norR, A_norR, O_norR] = deal(3, 50, 50);
[xNorR, yNorR] = tern(C_norR, A_norR, O_norR);
text(xNorR, yNorR, 'NORITE', 'FontSize', 9, 'FontWeight','bold', ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle','Rotation', 60);