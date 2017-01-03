
% wrapper to visualise LFP slowdown models

function [lag] = LFP_slowdown_wrapper(S,p)




% ci = 1;
% clear data
% lags = [70 95 113 168];
% for i = [1 8 9 10]
% 	data(ci).stimulus = PID(30e3:55e3,find(paradigm==i,1,'last'));
% 	data(ci).response = lags(ci)*ones(length(data(ci).stimulus),1);
% 	ci = ci+1;
% end


% list parameters
p.A;
p.B;
p.adap_tau;
p.K_tau;

% bounds
lb.A = 0;
lb.B = 0;
lb.adap_tau = 100;
lb.K_tau = 10;

R = asNL_euler(S,p);
x = S(5e3:end); x = x - mean(x); x = x/std(x);
y = R(5e3:end); y = y - mean(y); y = y/std(y);
d = finddelay(x,y);
if isnan(d) || isinf(d)
	d = 0;
end
lag = ones(length(S),1)*d;
