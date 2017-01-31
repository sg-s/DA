% non-adapting, non-parametric NL model

function [R,K,a] = NL(S,p)

% list parameters for clarity 
p.k_D;
p.n; 

% bounds
lb.k_D = 1e-3;
lb.n = 4;
ub.n = 4;

% decouple stimulus
T = S(:,2);
S = S(:,1);

% pass stimulus through nonlinearity 
a = 1./(1 + (p.k_D./S).^p.n);

% back out filter non-parametrically 
K = fitFilter2Data(a,T,'filter_length',1e3,'offset',200);
K = K(100:end-100);
filtertime = (1:length(K)) - 100;
time = 1:length(S);

% normalise filter correctly 
fp = convolve(time,a,K,filtertime);
K = K/(nanstd(fp)/nanstd(S)); % normalise correctly 
fp = convolve(time,a,K,filtertime);


rm_this = isnan(fp) | isnan(T);


x = fp(~rm_this);
y = T(~rm_this); 

ff = fit(x(:),y(:),'poly1');
R = fp*ff.p1;



