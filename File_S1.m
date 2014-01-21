%% This is a script to show standard use of the DA_integrate routine.

%% This work is licensed under a Creative Commons Attribution-ShareAlike 3.0 Unported License.
%% CC-BY-SA
%% Damon A. Clark, 2013

%% First, this is the form of the parameter set

% there are eight parameters here
p.A = 1; % alpha in model
p.B = 0.2; % beta in the model
p.C = 0.5; % gamma in the model
p.tau_r = 30; % these are the timescales
p.tau_y = 50;
p.n_y = 2;
p.tau_z = 60;
p.n_z = 3;

%% Next, generate the stimulus to be used

ints = 2.^[2:8];
S = reshape([zeros(length(ints),200),ints'*ones(1,1000),zeros(length(ints),500)]',1,[]);
% this is an increasing set of pulses

%% Now, call the integration routine

R = DA_integrate(S,p);

%% Compare with tau_r set to 0
p.tau_r=0;
R0 = DA_integrate(S,p);

%% Plot them up

figure;
subplot(2,1,1);
plot(S,'k');
ylabel('stimulus amplitude');
subplot(2,1,2); hold on;
plot(R,'k');
plot(R0,'r');
legend('tau_r = 30','tau_r = 0');
ylabel('response amplitude');
xlabel('time (a.u.)');




