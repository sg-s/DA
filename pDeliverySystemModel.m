% DeliverySystemModel.m
% this is a parameterised version of model of the odour delivery system
%
% created by Srinivas Gorur-Shandilya at 1:41 , 26 February 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function PID_pred = pDeliverySystemModel(s,p)



% core parameters -- will not be manipulated
MFC_Scale = 100; % mL/min/V
TotalFlow = 2000; % mL/min

% specify parameters and bounds
lb.n_PID = 1;
lb.n_MFC = 1;
lb.hill2_n = 1;
lb.hill1_n = 1;
lb.tau_MFC = 1;
lb.tau_PID = 1;
lb.hill2_K = 0;
lb.hill1_K = 0;

% MFC filter
p.A_MFC; 
p.tau_MFC; % for filter
p.n_MFC;

p.hill1_A;
p.hill1_K;
p.hill1_n;

% PID filter
p.A_PID; 
p.tau_PID; % for filter
p.n_PID;

p.hill2_A;
p.hill2_K;
p.hill2_n;


p_PID.A = p.A_PID; 
p_PID.tau = p.tau_PID;
p_PID.n = p.n_PID;

p_MFC.A = p.A_MFC; 
p_MFC.tau = p.tau_MFC;
p_MFC.n = p.n_MFC;

% pass through MFC filter
K_MFC = filter_gamma(p_MFC);
s = filter(K_MFC,1,s);

% and MFC hill function 
p_hill1.A = p.hill1_A;
p_hill1.k = p.hill1_K;
p_hill1.n = p.hill1_n;
s = hill(s,p_hill1);

% convert this into a dilution
dil = s*MFC_Scale;
dil = dil./(TotalFlow+dil);

% pass through PID filter
K_PID = filter_gamma(p_PID);
PID_pred = filter(K_PID,1,dil);

% pass through PID nonlinearity 
p_hill2.A = p.hill2_A;
p_hill2.k = p.hill2_K;
p_hill2.n = p.hill2_n;
PID_pred = hill(PID_pred,p_hill2);








