% fully parametric adapting NLN model
% this version is more DA-like, and uses a mono-lobed gamma filter to change K_D
% 

function [R,Ky,Kz,k_D,x] = aNL2(S,p)

% list parameters for clarity 

% input nonlinearity parameters
p.k0;
p.tau_z; % timescale of k_D update
p.B; % prefactor on k_D update
p.n_z;
p.n; % steepness of the Hill function

% parameters for response filter 
p.tau_y;
p.n_y;

% parameters for output scale
p.C; 
p.D;

% bounds for response
lb.k0 = 0;
lb.A = 0;
lb.C = 0;
lb.tau_y = 5;
lb.n_y = 1;
lb.n = 1;

ub.tau_y = 100;
ub.n = 8; 
ub.n_y = 10;

% bounds for adaptation 
lb.B = 0;
lb.tau_z = 1;
lb.n_z = 1;

ub.B = Inf;
ub.tau_z = 1000;
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
Ky = generate_simple_filter(p.tau_y,p.n_y,t);
temp = abs(Ky); temp = temp/max(temp);
Ky = Ky(1:find(temp>1e-2,1,'last'));

% pass through filter
R = filter(Ky,1,x);

% output nonlinearity 
R = R+p.D;
R = R*p.C;


function f = generate_simple_filter(tau,n,t)
f = t.^n.*exp(-t/tau); % functional form in paper
f = f/tau^(n+1)/gamma(n+1); % normalize appropriately

