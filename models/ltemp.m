function [R,H1,H2] = ltemp(S,p)

p.theta;
p.n;
p.b4;
p.b5;


H1 = (p.b4*S.^p.n)./(S.^p.n + p.theta.^p.n);
H2 = p.b5*S;

R = H1 + H2;