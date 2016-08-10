%% modified version of DA model designed to slow down
% 

function [R,y,z] = LFPmodel2(S,p)



% list parameters for legibility
p.A;
p.B;
p.tau_A;
p.tau_y;
p.tau_z;
p.tau_r;
p.s0;

t = 0:3000; % filters to be this long; don't use n*tau longer than a few hundred ms in this case...
% Kz and Ky are the filters described in equations 12 and 13
Ky = generate_simple_filter(p.tau_y,2,t);
Kz = generate_simple_filter(p.tau_z,2,t);

S = S + p.s0;

% y and z are the stimulus convolved with the filters Ky and Kz
y = filter(Ky,1,S);
z = filter(Kz,1,S);

% in the case that tau_r is 0, don't have to integrate; this approximation can
% also be used when tau_r/(1+p.B*z) << other time scales in y or z, for all or most z
if p.tau_r == 0
    R = p.A*y./(1+p.B*z);
    R = R(:);
    return;
end

% set up some variables to pass to the integration routine
pass.y = y;
pass.z = z;
pass.A = p.A;
pass.B = p.B;
pass.tau_A = p.tau_A;
pass.tau_r = p.tau_r;


% these options seem to work well, but can be modified, of course
opts = odeset('reltol',1e-6,'abstol',1e-5,'MaxStep',25);
T0 = [1,length(z)];
X0 = 0; % start from 0

% the equations can be quite stiff, and this ode solver works quite well
[tout,xout]=ode23t(@dxdt,T0,X0,opts,pass);

% linearly interpolate the output to time intervals of the inputs
R = interp1(tout,xout,1:length(z),'linear');
R = R(:);




function f = generate_simple_filter(tau,n,t)
f = t.^n.*exp(-t/tau); % functional form in paper
f = f/tau^(n+1)/gamma(n+1); % normalize appropriately

function dx = dxdt(t,x,pass)

zt = interp1(1:length(pass.z),pass.z,t,'linear'); % this is really slow
yt = interp1(1:length(pass.y),pass.y,t,'linear');

tau_r = (1+pass.tau_A*zt)*(1+pass.B*zt)/(1+zt)*pass.tau_r;

% this is the DA model equation at time t
dx = 1/tau_r * (pass.A*yt - (1 + pass.B*zt) * x);




