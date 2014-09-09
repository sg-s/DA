
function [Rguess, K] = SolveNLNModel(x,stimulus)	

% unpack parameters of model
tau1 = x(4); n1 = x(5); tau2 = x(6); n2 = x(7);   % bi-lobed filter block

% pass stimulus through input non-linearity
a = hill(x(1:3),stimulus);

% pass output through filter
t = 1:300;
K = make_bilobe_filter(tau1,n1,tau2,n2,t);
a = filter(K,1,a);

% pass it through a rectifier
a(a<0) = 0;


% pass filtered output through output non-linearity
Rguess = hill(x(8:10),a(:));

% throw out negative values
Rguess(Rguess<0) = 0;



