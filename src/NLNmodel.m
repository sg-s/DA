% simple NLN model with a single parameter (K_D)
% since we find the filter non-parameterically, we have to give it the response too
% S is a Tx2 matrix, where S(:,2) is the response we want to fit it to
% 
function [R,K] = NLNmodel(S,p)

% list parameters for clarity 
p.k_D;
p.n;


% bounds
lb.n = 1;
lb.k_D = 1e-2;

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


R(R<0) = 0;

rm_this = isnan(R) | isnan(T);

ff = fit(R(~rm_this),T(~rm_this),'poly1');
R = R*ff.p1;


