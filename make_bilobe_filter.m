function f = make_bilobe_filter(tau1,n1,tau2,n2,t)
f1 = t.^n1.*exp(-t/tau1); % functional form in paper
f1 = f1/tau1^(n1+1)/gamma(n1+1); % normalize appropriately


f2 = t.^n2.*exp(-t/tau2); % functional form in paper
f2 = f2/tau2^(n2+1)/gamma(n2+1); % normalize appropriately

f = f1 - f2;
f = f/max(f);