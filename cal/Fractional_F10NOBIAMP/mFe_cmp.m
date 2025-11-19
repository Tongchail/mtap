%% prepare workspace
clear; close all;

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
%% Calculate the mfe components
  b = [32.6154    4.5000    1.5192   36.3654   10.8558    9.7500    0.5481         0    1.9231    1.9231]; % Pietruszka.,2023
% b = [33.9020    5.1765    4.7255   23.3039    5.1961   12.5392    1.7353    0.6863   10.7745    1.9608]; % Honour.,2019
% b = [30.4020    1.9902    4.9412   36.6569    2.6078    8.0294    0.7059    0.4902   12.2157    1.9608]; % Hou.,2018
% b = [32.3065	  3.8889	3.7286	 32.1087	6.2199	 10.1062	0.9964	  0.3922	8.3044	  1.9482]; % Average ref.

A = [ 0    0.2200    4.0000   89.0200    6.7600         0         0         0         0         0
      0    5.3800    1.7300   89.1600    3.7300         0         0         0         0         0
      0    1.0900    0.8400   94.9300    3.1400         0         0         0         0         0

52.9200         0    4.0000   12.2000   28.9300    1.9500         0         0         0         0
53.3200         0    1.8900   18.0300   25.3200    1.4400         0         0         0         0
55.0000         0    0.9000   15.3700   28.1400    0.5900         0         0         0         0

53.4200         0    1.7900    7.4700   16.4000   20.1400    0.7800         0         0         0
54.5600         0    0.6500    9.3900   15.2300   18.8100    1.3600         0         0         0
55.4400         0    0.3100    8.3500   14.7400   19.0800    2.0800         0         0         0

      0   51.0322         0   42.5365    6.4313         0         0         0         0         0

      0         0         0         0         0   56.0000         0         0   42.0000         2

      0         0         0         0         0         0         0         0         0  100.0000];

A = A';

C = lsqregcmp(A, b,[0,0,0,0.0])*100;  % The last two are irrelevant here and can be set to 0. The first number makes the distribution across all end members more even if it is increased (or try [0,0.25,0.5,0.75]). The second one minimises the jump in fitted fractions between adjacent end members.
mFe_oxdfit = C*A.'/100;
res_mFe = mFe_oxdfit-b;
resnorm = rms(res_mFe);

disp('The component of mfe is:')
disp(C)

%   mmt	    tmt	       mgt	     chy	  fhy	     hyp     	mau     fau	      aug      ilm	    apt        wat