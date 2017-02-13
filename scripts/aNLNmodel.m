% adaptive NLN model
% minimal model with non-parametric filter estimation
% just like NLNmodel, but here we allow for adaptation by changing k_D in a feed-forward way from S

function [R,K,k_D,Kz] = aNLNmodel(S,p)

% list parameters for clarity 
p.n;

% adaptation
p.k0;
p.B;
p.tau;
p.n_z;

% bounds
lb.n = .5;
lb.n_z = 0;
lb.tau = 1;
lb.k0 = 0;
lb.B = 0;


ub.n = 4;

T = S(:,2); % target response
S = S(:,1); 

% build the adaptive filter
t = 0:(length(S)/10);
Kz = generate_simple_filter(p.tau,p.n_z,t);
temp = abs(Kz); temp = temp/max(temp);
Kz = Kz(1:find(temp>1e-2,1,'last'));
Kz = Kz/max(Kz);

% generate the dynamically updating k_D 
k_D = p.k0 + p.B*filter(Kz,1,S);

% pass stimulus through input non-linearity 
x = (S.^p.n)./(S.^p.n+k_D.^p.n);

% fit a filter from x to T
K = fitFilter2Data(x,T,'reg',1,'filter_length',700,'offset',100);

K = K(50:end-50);
filtertime = 1e-3*(1:length(K)) - 50e-3;
time = 1e-3*(1:length(T));
R = convolve(time,x,K,filtertime);


R(R<0) = 0;


rm_this = isnan(R) | isnan(T);

ff = fit(R(~rm_this),T(~rm_this),'poly1');
R = R*ff.p1;



function f = generate_simple_filter(tau,n,t)
f = t.^n.*exp(-t/tau); % functional form in paper
f = f/tau^(n+1)/gamma(n+1); % normalize appropriately
