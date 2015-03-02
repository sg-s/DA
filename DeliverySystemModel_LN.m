% DeliverySystemModel_LN.m
% this is a model of the odour delivery system
% this model is constructed a little strangely because
% the parameters of the model have to be specified in the base workspace. 
% the reason for doing this is because this model is called by another function, that
% we are trying to numerically optimise, and don't want to optimise the parameters of this model
% while doing so
%
% created by Srinivas Gorur-Shandilya at 1:41 , 26 February 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function PID_pred = DeliverySystemModel_LN(s)

PID_pred = DeliverySystemModel(s);

% pass through nonlinearity 
p_LN = evalin('base', 'p_LN');
PID_pred = hill(PID_pred,p_LN);