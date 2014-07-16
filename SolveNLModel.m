% SolveNLModel.m
% solves a Nonlinear-Linear model, i.e, a linear filter that operates on a static input non-linearity. 
% Usage:
% 1) generate output given a filter and a nonlinearity:
% [Rguess,K] = SolveNLModel(x,stimulus,R);
% where x is a 3-element vector specifying the shape of the logistic input nonlinearity 
% stimulus is a long vector
% and R is a shorter vector that is the filter, in the same time units as stimulus. 
%
% if R is as long as stimulus, it is interpreted as the response, and a filter K is fit to the data. 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [Rguess, K] = SolveNLModel(x,stimulus,R)	
if ~nargin
	help SolveNLModel
	return
end

% pass stimulus through input non-linearity
a = logistic(x,stimulus);

if length(R) == length(stimulus)
	% fit filter.
	K = FitFilter2Data(a,R,[],'reg=1;');

else
	% use provided filter. 
	K = R;

end

Rguess = filter(K,1,a);

% rectify it
Rguess(Rguess<0)=0;



