% aNL.m
% adaptive LN model, where we have mean and contrast adaptation baked in
% 
% created by Srinivas Gorur-Shandilya at 1:14 , 16 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [R,Ky,Kz,k_D,x] = aNL(S,p)

% list parameters for clarity 

% input nonlinearity parameters
p.k0;
p.tau_z; % timescale of k_D update
p.B; % prefactor on k_D update
p.n_z;
p.n; % steepness of the Hill function

% parameters for response filter 
p.tau1;
p.tau2;
p.n_y;
p.A;

% parameters for output scale
p.C; 

% bounds for response
lb.k0 = 1e-6;
lb.A = 0;
lb.C = 0;
lb.tau1 = 5;
lb.tau2 = 1;
lb.n_y = 2;
lb.n = 1;

ub.A = 0;
ub.tau1 = 1e3;
ub.tau2 = 1;
ub.n = 1; 
ub.n_y = 2;

% bounds for adaptation 
lb.B = 0;
lb.tau_z = 2e3;
lb.n_z = 2;

ub.B = 1e6;
ub.tau_z = 1e4;
ub.n_z = 2;




% use filters that are long enough, just about
t = 0:(length(S)/10);
Kz = generate_simple_filter(p.tau_z,p.n_z,t);
temp = abs(Kz); temp = temp/max(temp);
Kz = Kz(1:find(temp>1e-2,1,'last'));

% generate the dynamically updating k_D 
k_D = p.k0 + p.B*filter(Kz,1,S);

% pass stimulus through input non-linearity 
x = (S.^p.n)./(S.^p.n+k_D.^p.n);

% make the parametric  respone filter
p.n = p.n_y;
if p.A == 0
	p.tau = p.tau1;
	p.A = 1;
	Ky = filter_gamma(t,p);
else
	Ky = filter_gamma2(t,p);
end
temp = abs(Ky); temp = temp/max(temp);
Ky = Ky(1:find(temp>1e-2,1,'last'));

% pass through filter
R = filter(Ky,1,x);

% output nonlinearity 
R = R*p.C;


function f = generate_simple_filter(tau,n,t)
f = t.^n.*exp(-t/tau); % functional form in paper
f = f/tau^(n+1)/gamma(n+1); % normalize appropriately

