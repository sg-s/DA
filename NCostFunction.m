% NCostFunction.m
% computes the cost function for the Nonlinear-Linear model.
% meant to be called by FitNLModel.m
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function cost = NCostFunction(x,data)

Rguess = hill(x,data.stimulus);

% calculate the cost
cost = Cost2(data.response(300:end),Rguess(300:end));
