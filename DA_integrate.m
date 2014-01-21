function R = DA_integrate(S,p)

%% function takes argument of stimulus and a parameter set and implements the DA
%% model as described in Clark et al., 2013
%% the parameter set is as described in DA_model_script.m

%% This work is licensed under a Creative Commons Attribution-ShareAlike 3.0 Unported License.
%% CC-BY-SA
%% Damon A. Clark, 2013

t = [0:3000]; % filters to be this long; don't use n*tau longer than a few hundred ms in this case...
% Kz and Ky are the filters described in equations 12 and 13
Ky = generate_simple_filter(p.tau_y,p.n_y,t);
Kz = p.C*Ky + (1-p.C) * generate_simple_filter(p.tau_z,p.n_z,t);

% y and z are the stimulus convolved with the filters Ky and Kz
y = filter(Ky,1,S);
z = filter(Kz,1,S);

% in the case that tau_r is 0, don't have to integrate; this approximation can
% also be used when tau_r/(1+p.B*z) << other time scales in y or z, for all or most z
if p.tau_r == 0
    R = p.A*y./(1+p.B*z);
    return;
end

% set up some variables to pass to the integration routine
pass.y = y;
pass.z = z;
pass.tau_r = p.tau_r;
pass.A = p.A;
pass.B = p.B;

% these options seem to work well, but can be modified, of course
opts=odeset('reltol',1e-6,'abstol',1e-5,'MaxStep',25);
T0 = [1,length(z)];
X0 = 0; % start from 0
% X0 = p.A*y(1)/(1+p.B*z(1)); % start from steady state of first inputs

% the equations can be quite stiff, and this ode solver works quite well
[tout,xout]=ode15s(@dxdt,T0,X0,opts,pass);

% linearly interpolate the output to time intervals of the inputs
R = interp1(tout,xout,[1:length(z)],'linear');


% uncomment this if you want to check out the filters
% figure; hold on;
% plot(cumsum(Ky)); 
% plot(cumsum(Kz));





function f = generate_simple_filter(tau,n,t)

f = t.^n.*exp(-t/tau); % functional form in paper
f = f/tau^(n+1)/gamma(n+1); % normalize appropriately

function dx = dxdt(t,x,pass)

B = pass.B;
A = pass.A;
tau_r = pass.tau_r;
% can also use qinterp1 here for some speed up; see mathworks website
% maybe he means interp1q??? -- can't find qinterp1
% the reason this is happening is because ode15s calls this function and wants to know the value
% of z(t) and y(t) at arbitrary points, but we only know it at integer values of t (in ms)
zt = interp1([1:length(pass.z)],pass.z,t,'linear'); % this is really slow
yt = interp1([1:length(pass.y)],pass.y,t,'linear');
%zt = interp1q([1:length(pass.z)],pass.z,t); -- gives an error
%yt = interp1q([1:length(pass.y)],pass.y,t);

% this is the DA model equation at time t
dx = 1/tau_r * (A*yt - (1 + B*zt) * x);




