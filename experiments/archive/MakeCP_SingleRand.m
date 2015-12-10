m = 0:.2:2;
r = 2; % Volts, range

ControlParadigm.Name = [];
ControlParadigm.Outputs = [];

dt  =1e-4;
T = 60;
t_on = 2;
t_off = 58;
rand_on = 12;
rand_off = 55;

tc = .05; % 50ms

nsteps= (rand_off-rand_on)/tc;
noise = rand(nsteps,1);

for i = 1:length(m)
    
    n = strcat('mean-',mat2str(m(i)));
    ControlParadigm(i).Name = n;
    
    ControlParadigm(i).Outputs= zeros(2,T/dt);
    
    % turn the steady state on
    ControlParadigm(i).Outputs(1,t_on/dt:t_off/dt) = m(i);
    
    % set to a random number every tc
    rmin = m(i) - r;
    rmax = r + m(i);
    
    a = rand_on/dt;
    z = a + tc/dt;
    for j = 1:nsteps
        ControlParadigm(i).Outputs(1,a:z) = noise(j)*r + m(i);
        
        % increment
        a = z;
        z = a + tc/dt;
    end
    
    % clip
    ControlParadigm(i).Outputs(1,ControlParadigm(i).Outputs(1,:) > 5) = ControlParadigm(i).Outputs(1,ControlParadigm(i).Outputs(1,:) > 5) - 5;
    ControlParadigm(i).Outputs(1,ControlParadigm(i).Outputs(1,:) < 0) = ControlParadigm(i).Outputs(1,ControlParadigm(i).Outputs(1,:) < 0) + 5;
    
    % add the main air
    ControlParadigm(i).Outputs(2,:) = 1;
    ControlParadigm(i).Outputs(2,1) = 0;
    ControlParadigm(i).Outputs(2,end) = 0;
    
    

    
end