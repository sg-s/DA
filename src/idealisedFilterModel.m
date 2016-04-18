function [R,K] = idealisedFilterModel(S,p)

% list parameters for clarity
p.tau1;
p.tau2;
p.s1;
p.s2;

p.offset;

% bounds
lb.tau1 = 1;
lb.tau2 = 2;
lb.s1 = 0;
lb.s2 = 0;

% construct the filter
K = zeros(1e3,1);
t1 = round(p.tau1);
t2 = round(p.tau2);

K(1:t1) = linspace(0,-p.s1,t1);
K(t1+1:t2) = linspace(-p.s1,p.s2,t2-t1);
overflow = ceil((t2-t1)/2);
K(t2+1:t2+overflow) = linspace(p.s2,0,overflow);

R = filter(K,1,S);
R(isnan(R)) = 0;

R = R + p.offset;