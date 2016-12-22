% fully parametric adapting NLN model
% this version is more DA-like, and uses a mono-lobed gamma filter to change K_D
% 

function [R,Ky,Kz,k_D] = aNLN2(S,p)

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
lb.tau2 = 5;
lb.n_y = 2;
lb.n = 1;

ub.tau1 = 100;
ub.tau2 = 200;
ub.n = 32; 
ub.n_y = 2;

% bounds for adaptation 
lb.B = 0;
lb.tau_z = 1;
lb.n_z = 1;

ub.B = Inf;
ub.tau_z = 1e3;
ub.n_z = 2;




% see https://github.com/sg-s/DA/issues/114 for an explanation of the following
filter_length = 4*max([p.n_z*p.tau_z  p.n_y*p.tau1 p.n_y*p.tau2]);
if filter_length < length(S)/10
else
	filter_length = length(S)/10; % ridiculously long filters
end
t = 0:filter_length; 


% generate the dynamically updating k_D 
Kz = generate_simple_filter(p.tau_z,p.n_z,t);
k_D = p.k0 + p.B*filter(Kz,1,S);

% pass stimulus through input non-linearity 
x = (S.^p.n)./(S.^p.n+k_D.^p.n);

% make the parametric  respone filter
p.n = p.n_y;
Ky = filter_gamma2(t,p);

% pass through filter
R = filter(Ky,1,x);

% output nonlinearity 
R(R<0) = 0;
R = R*p.C;


function f = generate_simple_filter(tau,n,t)
f = t.^n.*exp(-t/tau); % functional form in paper
f = f/tau^(n+1)/gamma(n+1); % normalize appropriately

