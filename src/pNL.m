% fully parametric non-adapting NL model
% 
% since we find the filter non-parameterically, we have to give it the response too
% S is a Tx2 matrix, where S(:,2) is the response we want to fit it to
% 
function [R] = pNL(S,p)

% list parameters for clarity 

% input nonlinearity parameters
p.Hill_n;
p.Hill_K;

% parameters for filter 
p.tau1;
p.tau2;
p.n;
p.A;

% parameters for output scale
p.C; 

% lower bounds
lb.A = 0;
lb.B = 0;
lb.C = 0;
lb.tau1 = 5;
lb.tau2 = 20;
lb.n = 1; 
lb.Hill_n = 8;
lb.Hill_K = 1e-3;

% upper bounds
ub.tau1 = 200;
ub.tau2 = 400;
ub.n = 10;
ub.Hill_n = 32; 


% pass stimulus through input non-linearity 
x = (S.^p.Hill_n)./(S.^p.Hill_n+p.Hill_K.^p.Hill_n);

% make the parametric filter
K = filter_gamma2(1:1e3,p);

% pass through filter
R = filter(K,1,x);

% output nonlinearity 
R(R<0) = 0;
R = R*p.C;

