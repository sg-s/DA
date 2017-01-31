
% wrapper to visualise LFP slowdown models

function [lag] = LFP_slowdown_wrapper(S,p)


% nsteps = 5;
% b = logspace(log10(.2),log10(2),nsteps);
% S = repmat([ones(10e3,1); 2*ones(10e3,1); ones(5e3,1)],1,nsteps);
% for i = 1:nsteps
% 	S(:,i) = S(:,i)*b(i);
% end

% ci = 1;
% clear data
% for i = [1:max(paradigm)]
% 	data(ci).stimulus = PID(30e3:50e3,find(paradigm==i,1,'last'));
% 	l = LFP_lags(paradigm==i); l(l<0) = NaN;
% 	data(ci).response = nanmean(l)*ones(length(data(ci).stimulus),1);
% 	ci = ci+1;
% end
% data(isnan(mean([data.response]))) = [];


% list parameters
p.E0;
p.E1;
p.w_minus;
p.w_plus;

p.tau_adap;
p.kT; 


% output filter
p.K_tau;

% bounds 
lb.E0 = 0;
lb.E1 = 0;
lb.w_minus = 0;
lb.w_plus = 0;
lb.tau_adap = .5;
lb.kT = 0;
lb.K_tau = 0;

R = asNL3(S,p);
x = S(5e3:end); x = x - mean(x); x = x/std(x);
y = R(5e3:end); y = y - mean(y); y = y/std(y);
d = finddelay(x,y);
if isnan(d) || isinf(d)
	d = 0;
end
lag = ones(length(S),1)*d;
