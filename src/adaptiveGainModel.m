% adaptiveGainModel.m
% model that attempts to fit an adaptive gain term (the denom. of the DA model) to data
% 
% created by Srinivas Gorur-Shandilya at 2:09 , 03 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [R,Kg] = adaptiveGainModel(S,p)

switch nargin
case 0
	help adaptiveGainModel
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

% declare parameters for clarity
p.A;
p.B;
p.tau;
p.n;

% bounds
lb.A = eps;
lb.B = 100;
lb.tau = 1;
lb.n = 1;
ub.tau = 100;
ub.n = 10;

t = 0:(length(S)/10); 
Kg = generate_simple_filter(p.tau,p.n,t);
% clip filter to interesting bits
[m,loc]=max(Kg);
z = loc+find(Kg(loc:end)<m/100,1,'first'); % accept a 1% error
if ~isempty(z)
	Kg = Kg(1:z);
end

fp = S(:,1);
S = S(:,2);

g = filter(Kg,1,S);
R = (p.A*fp)./(1+p.B.*g);


function f = generate_simple_filter(tau,n,t)
f = t.^n.*exp(-t/tau); % functional form in paper
f = f/tau^(n+1)/gamma(n+1); % normalize appropriately