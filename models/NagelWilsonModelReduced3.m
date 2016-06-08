% NagelWilsonModelReduced3.m
% reduced version of Nagel and Wilson's model where we assume that binding timescales are much shorter than conformational change timescales
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function dy = NagelWilsonModelReduced3(t,y,time,odor,p)

% calculate the odor at the time point
O = interp1(time,odor,t); % Interpolate the data set (ft,f) at time t

% receptor-ligand binding -- this has been reduced to one equation
R = y(1);
C = y(2);
D = y(3);

% now calculate Rx and ORx
Rx = (1 - R - O*R*p.kb)/(1 + p.theta*p.kb*O);
ORx = Rx*p.theta*p.kb*O;
assert(Rx>0,'Rx is -ve')
assert(ORx>0,'ORx is -ve')


% solve the odes
dR = -p.ka*p.sa*R + p.sa*Rx;
dC = p.ko*(1/D)*(1-C)*(Rx+ORx) - p.kc*C;
dD = p.A*C - p.B*D;

dy = [dR; dC; dD];


