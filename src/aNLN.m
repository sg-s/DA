% fully parametric adapting NLN model
% 
% since we find the filter non-parameterically, we have to give it the response too
% S is a Tx2 matrix, where S(:,2) is the response we want to fit it to
% 
function [R,K,k_D] = aNLN(S,p)

% list parameters for clarity 

% input nonlinearity parameters
p.k0;
p.tau; % timescale of k_D update
p.B; % prefactor on k_D update

% parameters for filter 
p.tau1;
p.tau2;
p.n;
p.A;

% parameters for output scale
p.C; 

% hard bounds
lb.k0 = 1e-3;
lb.n = 1;
lb.tau = 1;
lb.A = 0;
lb.B = 0;
lb.C = 0;
lb.tau1 = 10;
lb.tau2 = 20;

% upper bounds
ub.tau = 1e3;
ub.tau2 = 200;
ub.n = 4; % anything more is unreasonable 


% generate the dynamically updating k_D 
K = ones(ceil(p.tau),1);
k_D = p.k0 + p.B*filter(K,length(K),S);

% pass stimulus through input non-linearity 
x = (S.^p.n)./(S.^p.n+k_D.^p.n);

% make the parametric filter 
K = filter_gamma2(1:1e3,p);

% pass through filter
R = filter(K,1,x);

% output nonlinearity 
R(R<0) = 0;
R = R*p.C;
