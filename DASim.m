% DASim.m
% simulates synthetic stimuli, using a exponentiated gaussian.
% then simulates output using a DA Model.
% then performs gain analysis. 
% this is meant to work with Manipulate.m
% calling DASim with no arguments will make the figure, and return handles to the relevant axes
% calling DASim with an argument assumes you have the figure ready, and are willing to accept data to plot
% 
% created by Srinivas Gorur-Shandilya at 6:28 , 06 February 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [] = DASim(p)

% check if the figure has been created
DASim_figure = [];
try DASim_figure = evalin('base', 'DASim_figure');
	DASim_plot = evalin('base', 'DASim_plot');
catch
	DASim_figure = figure('Units','pixels','Position',[10 100 1140 700],'Name','DASim','IntegerHandle','off','CloseRequestFcn',@closess,'Color','w'); hold on
	clf;
	DASim_plot(1) = axes('Units','pixels','Position',[35 560 700 122.5]);
	ylabel('Stimulus')
	DASim_plot(2) = axes('Units','pixels','Position',[35 350 700 192.5]);
	ylabel('Response')
	DASim_plot(3) = axes('Units','pixels','Position',[787.5 367.5 297.5 297.5]);
	title('DA Filters')
	DASim_plot(4) = axes('Units','pixels','Position',[35 35 315 280]);
	title('Filter')
	DASim_plot(5) = axes('Units','pixels','Position',[402.5 35 315 280]);
	xlabel('Linear Prediction')
	ylabel('Response')
	DASim_plot(6) = axes('Units','pixels','Position',[787.5 52.5 297.5 262.5]);
	xlabel('History Length')
	ylabel('Gain')

	% save to the base workspace
	assignin('base','DASim_figure',DASim_figure)
	assignin('base','DASim_plot',DASim_plot)
end

% some hard coded parameters
T = 60;
dt = 3e-3;
time = dt*(1:floor(T/dt));

% set up frozen noise
stream = RandStream.getGlobalStream;
savedState = [];
load('rand_state.mat');

stream.State = savedState;
noise = randn(1,ceil(T/dt));

% make the stimulus 
stimulus = p.baseline - 1 + exp(p.sigma*randn(floor(T/dt),1));
stimulus(stimulus<p.baseline) = p.baseline;

% filter by autocorrelation time
p.act = 1+round(p.act/dt);
stimulus = filtfilt(ones(p.act,1)/length(ones(p.act,1)),1,stimulus);
plot(DASim_plot(1),time,stimulus,'k')
ylabel(DASim_plot(1),'Stimulus')


% make the DA Model
p.tau_z = p.tau_z/dt;
p.tau_y = p.tau_y/dt; 
[R,~,~,Ky,Kz] = DAModelv2(stimulus,p); 
filtertime = dt*(1:length(Ky));
plot(DASim_plot(2),time,R,'k');
ylabel(DASim_plot(2),'Response')
plot(DASim_plot(3),filtertime,Ky,'g');
hold(DASim_plot(3),'on') 
plot(DASim_plot(3),filtertime,Kz,'m');
legend(DASim_plot(3),'K_y','K_z')
hold(DASim_plot(3),'off') 


% do the computationally difficult things only when the user is not actively doing things
event = [];
eval(strcat('event =','p','.event;')); % this bizzare construction is so that getModelParameters doesn't see this code. Yes, we're obfuscating code to hide from other code that was designed to read code. 
compute = 0;
if isempty(event)
	compute = 1;
elseif strcmp(event.EventName,'ContinuousValueChange')
	compute = 0;
elseif strcmp(event.EventName,'Action')
	compute = 1;
end

if compute
	% extract a linear filter
	[K, ~, filtertime] = FindBestFilter(stimulus,R,[],'regmax=1;','regmin=1;','filter_length=500;','offset=50;');
	filtertime = filtertime*dt;
	plot(DASim_plot(4),filtertime,K,'r')
	title(DASim_plot(4),'Filter')

	% make a linear prediction
	fp = convolve(time,stimulus,K,filtertime) + mean(R);
	hold(DASim_plot(2),'on')
	plot(DASim_plot(2),time,fp,'r');
	hold(DASim_plot(2),'off')

	% scatter plot
	ss = 3;
	plot(DASim_plot(5),fp(300:ss:end),R(300:ss:end),'.')
	xlabel(DASim_plot(5),'Linear Prediction')
	ylabel(DASim_plot(5),'Response')
end

% just so that getModelParameters can read it
p.act = 1;
p.n_y = 2; 
p.n_z = 2; p.A = 1; p.B = 1; p.C = 0; 
p.tau_z = 1; p.tau_y = 1;
p.s0 = 0;


% specify bounds for FitModel2Data
lb.A = 1; lb.B = 1; lb.C = 0 ; 
lb.tau_y = .01; lb.tau_z = .01;

ub.C = 1; ub.act = 0;

% extra bounds
lb.n_y = 2; lb.n_z = 2;
ub.n_y = 2; ub.n_z = 2;
lb.s0 = -5; ub.s0 = 1;
ub.tau_z = 1; ub.tau_y = 1;


function [] = closess(~,~)
	delete(DASim_figure)
	evalin('base','clear DASim_figure')
	evalin('base','clear DASim_plot')

end


end


