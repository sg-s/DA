% contrastLogisticModel
% the stimulus has two vectors:
% the first is the linear prediction
% and the second is the derivative of the stimulus
% 
% created by Srinivas Gorur-Shandilya at 1:49 , 16 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [R,Kg,k] = contrastLogisticModel(S,p)

switch nargin
case 0
	help contrastLogisticModel
	return
case 1
	error('Not enough input arguments')
case 2
	if size(S,2) < size(S,2)
		S = S';
	end
	assert(size(S,2) == 2,'Stimulus should have two columns')
	assert(isstruct(p),'Second argument should be a structure')
end


% for logistic function
p.x0;
p.A;

% for changing k:
p.k0;
p.n;
p.tau;
p.B;

% bounds
lb.tau = 1;
lb.n = 1;
lb.B = 0;
lb.A = 10;
lb.k0 = 0;

ub.tau = 400;
ub.n = 10;
ub.B = 1e6;
ub.x0 = 0;

filter_length = 4*(p.n*p.tau);
t = 0:filter_length; 
Kg = generate_simple_filter(p.tau,p.n,t);

LFP = S(:,1);
LFP_projected = S(:,2);

LFP(LFP<0) = 0;
Shat = filter(Kg,1,LFP);

k = p.k0./(1 + p.B*Shat);

R = logistic(LFP_projected,p.A,k,p.x0);


function f = generate_simple_filter(tau,n,t)
f = t.^n.*exp(-t/tau); % functional form in paper
f = f/tau^(n+1)/gamma(n+1); % normalize appropriately



