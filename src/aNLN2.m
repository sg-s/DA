% fully parametric adapting NLN model
% this version is more DA-like, and uses mono-lobed gamma filters everywhere
% 
% since we find the filter non-parameterically, we have to give it the response too
% S is a Tx2 matrix, where S(:,2) is the response we want to fit it to
% 
function [R,Ky,Kz,k_D] = aNLN2(S,p)

% list parameters for clarity 

% input nonlinearity parameters
p.k0;
p.tau_z; % timescale of k_D update
p.B; % prefactor on k_D update
p.n;
p.n_z;

% parameters for response filter 
p.tau_y;
p.n_y;
p.A;

% parameters for output scale
p.C; 

% lower bounds
lb.k0 = 1e-6;
lb.A = 0;
lb.B = 0;
lb.C = 0;
lb.tau_y = 5;
lb.tau_z = 20;
lb.n_y = 1;
lb.n_z = 5;

lb.n = 8; % constrained by the variance data?

% upper bounds
ub.tau_y = 100;
ub.tau_z = 200;
ub.n = 32; 
ub.n_z = 10;
ub.n_y = 10;

% see https://github.com/sg-s/DA/issues/114 for an explanation of the following
filter_length = 4*max([p.n_z*p.tau_z  p.n_y*p.tau_y]);
if filter_length < length(S)/10
else
	filter_length = length(S)/10; % ridiculously long filters
end
t = 0:filter_length; 


% generate the dynamically updating k_D 
Kz = generate_simple_filter(p.tau_z,p.n_z,t);
k_D = p.k0*0 + p.B*filter(Kz,1,S);

% pass stimulus through input non-linearity 
x = (S.^p.n)./(S.^p.n+k_D.^p.n);

% make the parametric filter
Ky = generate_simple_filter(p.tau_y,p.n_y,t);

% pass through filter
R = filter(Ky,1,x);

% output nonlinearity 
R(R<0) = 0;
R = R*p.C;


function f = generate_simple_filter(tau,n,t)
f = t.^n.*exp(-t/tau); % functional form in paper
f = f/tau^(n+1)/gamma(n+1); % normalize appropriately

