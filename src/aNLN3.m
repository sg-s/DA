% fully parametric adapting NLN model
% this version is more DA-like, and uses a mono-lobed gamma filter to change K_D
% this version has a minimum k_D that is attianable 

function [R,Ky,Kz,k_D,x] = aNLN3(S,p)

% list parameters for clarity 

% input nonlinearity parameters
p.tau_z; % timescale of k_D update
p.B; % strength of adaptation 
p.n_z;
p.n; % steepness of the Hill function
p.k_D_min;

% parameters for response filter 
p.tau1;
p.tau2;
p.n_y;
p.A;

% parameters for output scale
p.C; 

% bounds for response
lb.A = 0;
lb.C = 0;
lb.tau1 = 1;
lb.tau2 = 1;
lb.n_y = 1;
lb.n = 1;

ub.tau1 = 100;
ub.tau2 = 500;
ub.n = 4; 
ub.n_y = 2;

% bounds for adaptation 
lb.k_D_min = 1e-4;
lb.B = 0;
lb.tau_z = 100;
lb.n_z = 1;

ub.k_D_min = 10;
ub.B = 100;
ub.tau_z = 2e3;
ub.n_z = 2;




% use filters that are long enough, just about
t = 0:(length(S)/10);
Kz = generate_simple_filter(p.tau_z,p.n_z,t);
temp = abs(Kz); temp = temp/max(temp);
Kz = Kz(1:find(temp>1e-2,1,'last'));

% generate the dynamically updating k_D 
k_D = p.B*filter(Kz,1,S);
k_D(k_D<p.k_D_min) = p.k_D_min;

% pass stimulus through input non-linearity 
x = (S.^p.n)./(S.^p.n+k_D.^p.n);

% make the parametric  respone filter
p.n = p.n_y;
Ky = filter_gamma2(t,p);
temp = abs(Ky); temp = temp/max(temp);
Ky = Ky(1:find(temp>1e-2,1,'last'));

% pass through filter
R = filter(Ky,1,x);

% output nonlinearity 
R(R<0) = 0;
R = R*p.C;


function f = generate_simple_filter(tau,n,t)
f = t.^n.*exp(-t/tau); % functional form in paper
f = f/tau^(n+1)/gamma(n+1); % normalize appropriately

