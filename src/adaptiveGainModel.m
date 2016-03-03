% adaptiveGainModel.m
% model that attempts to fit an adaptive gain term (the denom. of the DA model) to data
% 
% created by Srinivas Gorur-Shandilya at 2:09 , 03 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function R = adaptiveGainModel(S,p)

switch nargin
case 0
	help adaptiveGainModel
	return
case 1
	error('Not enough input arguments')
case 2
	assert(isvector(S),'First argument should be a vector')
	assert(isstruct(p),'Second argument should be a structure')
end

% declare parameters for clarity
p.A;
p.B;
p.tau;
p.n;

% bounds
lb.A = 0;
lb.B = 0;
lb.tau = 0;
lb.n = 1;
ub.tau = 600;
ub.n = 10;

filter_length = 4*max(p.n*p.tau);
if filter_length < length(S)/10
else
	filter_length = length(S)/10; % ridiculously long filters
end
t = 0:filter_length; 

Kg = generate_simple_filter(p.tau,p.n,t);
g = filter(Kg,1,S);
R = (p.A*S)./(1+p.B*g);

function f = generate_simple_filter(tau,n,t)
f = t.^n.*exp(-t/tau); % functional form in paper
f = f/tau^(n+1)/gamma(n+1); % normalize appropriately