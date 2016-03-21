% 
% 
% created by Srinivas Gorur-Shandilya at 12:08 , 21 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function R = NLModel(S,p)

% parameters
p.K;
p.n;

lb.K = 0;
lb.n = .1;

R = S(:,1);
S = S(:,2);


S = real((S.^p.n)./(p.K.^p.n + S.^p.n));

K = fitFilter2Data(S,R,'offset',200,'filter_length',800);

K = K(100:end-100);
ft = (1:length(K)) - 100;
X = convolve(1:length(S),S,K,ft);

rm_this = isnan(X) | isnan(R);

ff = fit(X(~rm_this),R(~rm_this),'poly1');

R = ff(X);