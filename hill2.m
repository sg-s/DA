function r = hill2(x,xdata)
A = x(1);
K = x(2);
n = x(3);
r=real(A*(1./(1+(K./xdata).^n))) ;