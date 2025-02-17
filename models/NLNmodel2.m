% just like NLNmodel, but with an extra parameter that specifies the cutoff. 
% this also regresses the output to the response, fitting both the slope and any offset
function [R,K] = NLNmodel2(S,p)

% list parameters for clarity 
p.k_D;
p.n;

p.cutoff; 

% bounds
lb.n = .5;
lb.k_D = 0;
ub.n = 4;

T = S(:,2); % target response
S = S(:,1); 

% pass stimulus through input non-linearity 
x = (S.^p.n)./(S.^p.n+p.k_D^p.n);

% fit a filter from x to T
K = fitFilter2Data(x,T,'reg',1,'filter_length',700,'offset',100);

K = K(50:end-50);
filtertime = 1e-3*(1:length(K)) - 50e-3;
time = 1e-3*(1:length(T));
R = convolve(time,x,K,filtertime);

rm_this = isnan(R) | isnan(T);
ff = fit(R(~rm_this),T(~rm_this),'poly1');
R = ff(R);

R(R<p.cutoff) = 0;


