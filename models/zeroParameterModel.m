% zeroParameterModel.m
% 
% created by Srinivas Gorur-Shandilya at 11:21 , 17 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [R] = zeroParameterModel(S,p)


switch nargin
case 0
	if ~nargout
		help zeroParameterModel
		return
	end
case 1
	error('Not enough input arguments')
case 2
	assert(isvector(S),'First argument should be a vector')
	assert(isstruct(p),'Second argument should be a structure')
end

% specify parameters for clarity
p.tau_r; % response time scale
p.tau_m; % mean gain control timescale
p.tau_s; % contrast gain control timescale

p.B_m;   % mean gain control amount
p.B_s; 	 % contrast gain control amount

p.n0;    % max. Hill steepness
p.A; 	 % Hill maximum
p.K; 	 % Hill K_D

p.s0; 	 % stimulus offset term

% bounds
lb.tau_r = 1;   lb.tau_m = 1;   lb.tau_s = 1;
ub.tau_r = 1e2; ub.tau_m = 1e2; ub.tau_s = 1e2;

lb.B_m = 0; lb.B_s = 0;   lb.n0 = 1; lb.A = 1;  lb.K = 0;
ub.n0 = 10; ub.B_s = 10; 						ub.K = 1;


S = (S + p.s0); 

% see https://github.com/sg-s/DA/issues/114 for an explanation of the following
filter_length = 4*max([2*p.tau_m  2*p.tau_r]);
if filter_length < length(S)/10
else
	filter_length = length(S)/10; % ridiculously long filters
end
t = 0:filter_length; 

% Kz and Ky are the filters described in equations 12 and 13
Kr = generate_simple_filter(p.tau_r,2,t);
Km = generate_simple_filter(p.tau_m,2,t);

% y and z are the stimulus convolved with the filters Ky and Kz
y = filter(Kr,1,S);
z = filter(Km,1,S);

X = (y./(1+p.B_m*z));

% normalise this 
X = X - min(X);
X = X/max(X);

% calculate the contrast-dependent n
tau_s = round(p.tau_s);
Kg = [ones(tau_s,1); -ones(tau_s,1)];
Shat = filter(Kg,length(Kg),S);
Shat = abs(Shat);

n = p.n0./(1 + p.B_s.*Shat);

% pass through the Hill function
R = X.^n;
R = R./(p.K.^n + R.^n);

R = R*p.A;




function f = generate_simple_filter(tau,n,t)
f = t.^n.*exp(-t/tau); % functional form in paper
f = f/tau^(n+1)/gamma(n+1); % normalize appropriately



