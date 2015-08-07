% MakeCP_ShiftedGaussian.m
% make control paradigms for delivering a set of stimuli
% where the stimulus is drawn from a certain distribution with a certain
% standard deviation, and we move the mean of the distribution in each
% paradigm, keeping the standard deviation the same

rmin = -100; % defined relative to the mean flow rate
rmax = 100;
m = 0:20:100; % mL/min, flow through bbackground odor bottle
tc = .05; % switching time

% some fixed parameters
MaxFlows = [1000 500 500 200]; % for the MFCs
MaxValveFlow = 200; % mL/min

% time
dt = 1e-4;
T = 60;
f_on = 5;
f_off = 55;
b_on = 1;
b_off = T-1;


nsteps= (f_off-f_on)/tc;
noise = rand(nsteps,1);

ControlParadigm(1).Name = '';
ControlParadigm(1).Outputs = [];

for i = 1:length(m)
    n = strcat('mean-',mat2str(m(i)));
    ControlParadigm(i).Name = n;
    
    
    % initialise
    ControlParadigm(i).Outputs = zeros(7,T/dt);
    
    % set it to the mean
    b = m(i)/MaxFlows(4)*5;
    b0 = (MaxValveFlow-m(i))/MaxFlows(2)*5;
    ControlParadigm(i).Outputs(4,b_on/dt:b_off/dt) = b;
    ControlParadigm(i).Outputs(2,b_on/dt:b_off/dt) = b;
    
    this_rmin = rmin + m(i);
    this_rmax = rmax + m(i);
    
    if this_rmin < 0
        this_rmin = 0.1;
    end
    
    if this_rmax > MaxValveFlow
        this_rmax = MaxValveFlow;
    end
    lmin = (this_rmin);
    lmax = (this_rmax);
    
    a = f_on/dt;
    z = a + tc/dt;
    for j = 1:nsteps
        % pick a random number
        x = lmin+ noise(j)*(lmax-lmin);
        x = (x); % this is the flow rate we should set it to
        
        % covnert to a voltage and set it
        ControlParadigm(i).Outputs(4,a:z) = x/MaxFlows(4)*5;
        ControlParadigm(i).Outputs(2,a:z) = (MaxValveFlow - x)/MaxFlows(2)*5;
        
        % increment a,z
        a = a+tc/dt;
        z = z+tc/dt;
        
        
    end
    
    
end
