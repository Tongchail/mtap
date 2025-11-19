%% prepare workspace
clear all; close all;

addpath(genpath('../'));
addpath('../../cal')
addpath('../../src')
addpath('../../../unmix')
addpath('../../../unmix/src')
load ocean
Fs = {'FontSize',12};
FS = {'FontSize',15};
FL = {'FontSize',18};
Ms = {'MarkerSize',8};
MS = {'MarkerSize',10};
ML = {'MarkerSize',12};
TX = {'Interpreter','latex'};
TL = {'TickLabelInterpreter','latex'};
LB = {'Location','best'};
LO = {'Location','bestoutside'};

%% different calibration data
  cal_MtAp_750_new      % Pietruszka.,2023
 %cal_MtAp_750_honour   % Honour.,2019
 %cal_MtAp_750_hou      % Hou.,2018
 %cal_MtAp_750_avg      % Average ref.

%%  main 
b = [57.50    0.86    16.31   6.73   4.43    6.78    3.45    1.95    0.22    1.87]; % lowe SiO2 host andesite

f = b(end-1)./cal.cmp_oxd(6,end-1) * 0.5;  %  b(end-1)./cal.cmp_oxd(6,end-1)  (P₂O₅)
ba = (b - f.*cal.cmp_oxd(6,:))./(1-f);

A = cal.cmp_oxd;
A = A';

C = lsqregcmp(A, ba,[12,0,0,0.0])*100;
C(6) = f*100;
C = C./sum(C)*100;

start_oxdfit = C*A.'/100;
res = start_oxdfit-b;
resnorm = rms(res);

disp('The start components:')
disp(C)