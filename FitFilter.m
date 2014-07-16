% x = FitFilter(K)
% fits a bilobed filter made from two gamma distributions to an actual numerically computed filter 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function x = FitFilter(K)

if ~nargin
	help FitFilter
	return
end

t = 1:length(K);

x0 = [0.1 2 0.2 2];
lb = [0 0 0 0];
ub = [length(K) 100 length(K) 100];

psoptions = psoptimset('UseParallel',true,'CompletePoll', 'on', 'Vectorized', 'off','Display','iter','MaxIter',200,'MaxFunEvals',10000);
x = patternsearch(@(x) FFCostFunction(x,K,t),x0,[],[],[],[],lb,ub,psoptions);	

function [cost] = FFCostFunction(x,K,t)

% construct the filter
Kfit = make_bilobe_filter(x(1),x(2),x(3),x(4),t);

% calculate the cost
keyboard
cost = Cost2(Kfit,K)

