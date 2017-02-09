% divisiveGainModel.m
% this model simply handles the divisive gain term,
% and returns the expected gain at each mean stimulus
% you cannot use this to directly compute the response
% 
% created by Srinivas Gorur-Shandilya at 10:58 , 22 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [G,g_t] = divisiveGainModel(S,p)

% list parameters for legibility
p.n;
p.tau;
p.B;
p.A;

lb.n = 1;
lb.tau = 1;
lb.B = 0;

ub.n = 1;
ub.tau = 100;

% see https://github.com/sg-s/DA/issues/114 for an explanation of the following
filter_length = 4*(p.n*p.tau);
if filter_length < length(S)/10
else
	filter_length = length(S)/10; % ridiculously long filters
end
t = 0:filter_length; 

Kg = generate_simple_filter(p.tau,p.n,t);

temp = S; w = round(width(temp)/2);
S = temp(:,1:w);
fp = temp(:,w+1:end);

g_t = 0*S;

G = NaN(width(S),1);

for i = 1:width(S)
	g = 1./(1+ p.B*filter(Kg,1,S(:,i)));
	g_t(:,i) = g;
	y = fp(:,i).*g;
	% find the gain
	x = fp(:,i);
	rm_this = isnan(x) | isnan(y);
	rm_this(1:1e3) = true;
	x(rm_this) = []; y(rm_this) = [];
	ff = fit(x,y,'poly1');
	G(i) = ff.p1;
end

G = p.A*G;




function f = generate_simple_filter(tau,n,t)
f = t.^n.*exp(-t/tau); % functional form in paper
f = f/tau^(n+1)/gamma(n+1); % normalize appropriately



