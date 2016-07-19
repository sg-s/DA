%% implements the first order model:
% 0 = K0 x s(t) + (1 + K1 x s(t)) r(t)
% where x is a convolution
% to data specified by s(t) and r(t)
% and K0 and K1 are filters that are fit non-parametrically 
% this is basically a non-parametreric, simplified DA model

function [R] = firstOrderModel(S,filterset)

K0 = filterset.K0;
K1 = filterset.K1;
S = circshift(S,min(filterset.time)-1);	
R = filter(K0,1,S)./(1  + filter(K1,1,S));
