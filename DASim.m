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
marker_size=10;
marker_size2=24;
font_size=20;
plot_start = find(time>10,1,'first');

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
set(DASim_plot(1),'Xlim',[10 60])
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
set(DASim_plot(2),'Xlim',[10 60])


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
	fp = convolve(time,stimulus,K,filtertime);
	f=fit(fp(~(isnan(fp) | isnan(R))),R(~(isnan(fp) | isnan(R))),'poly1');
	fp = fp*f.p1;
	fp = fp+f.p2;

	hold(DASim_plot(2),'on')
	l=plot(DASim_plot(2),time,fp,'r');
	t2 = floor(length(fp)/2);
	r2 = strcat('r^2=',oval(rsquare(fp(t2:end),R(t2:end))));
	legend(l,r2)
	hold(DASim_plot(2),'off')
	set(DASim_plot(2),'Xlim',[10 60])

	% scatter plot
	ss = 3;
	plot(DASim_plot(5),fp(plot_start:ss:end),R(plot_start:ss:end),'.','MarkerSize',marker_size,'MarkerFaceColor',[0.9 0.9 0.9],'MarkerEdgeColor',[0.9 0.9 0.9])
	xlabel(DASim_plot(5),'Linear Prediction')
	ylabel(DASim_plot(5),'Response')
	hold(DASim_plot(5),'on')

	% compute shat(the smoothed stimulus)
	p.hl = floor(p.hl/dt)+1;
	shat = ComputeSmoothedStimulus(stimulus,round(p.hl));

	n = floor(sum(~isnan(R))*p.frac);


	% find times when smoothed stim is lowest x%
	this_shat = shat;
	this_shat(1:p.hl) = Inf; % the initial segment where we can't estimate shat is excluded
	this_shat(isnan(this_shat)) = Inf;
	this_shat(isnan(R)) = Inf;
	[~, t_low] = sort(this_shat,'ascend');
	t_low = t_low(1:n); % this is an index
	f_low = R(t_low);
	fp_low = fp(t_low);
	s_low = this_shat(t_low);
	t_low = time(t_low); % t_low is now a time. 
	 

	% find times when smoothed stim is highest x%
	this_shat = shat;
	this_shat(1:p.hl) = -Inf;
	this_shat(isinf(this_shat)) = -Inf;
	this_shat(isnan(R)) = -Inf;
	[~, t_high] = sort(this_shat,'descend');
	t_high = t_high(1:n);
	f_high = R(t_high);
	fp_high = fp(t_high);
	s_high = this_shat(t_high);
	t_high  = time(t_high);

	% scatter with colours
	clear l L
	l(1) = plot(DASim_plot(5),fp_low(1:ss:end),f_low(1:ss:end),'.','MarkerSize',marker_size,'MarkerFaceColor',[0.5 1 0.5],'MarkerEdgeColor',[0.5 1 0.5]);
	l(2) = plot(DASim_plot(5),fp_high(1:ss:end),f_high(1:ss:end),'.','MarkerSize',marker_size,'MarkerFaceColor',[1 0.5 0.5],'MarkerEdgeColor',[1 0.5 0.5]);
	hold(DASim_plot(5),'off')

	% remove NaNs
	censor_these = find(isnan(fp_high) + isnan(fp_low) + isnan(f_low) + isnan(f_high));
	f_high(censor_these) = [];
	fp_high(censor_these) = [];
	f_low(censor_these) = [];
	fp_low(censor_these) = [];
	t_low(censor_these) = [];
	t_high(censor_these) = [];

	% fit lines
	f = fit(fp_low(1:ss:end),f_low(1:ss:end),'poly1');
	L{1} = strcat('m=',oval(f.p1));
	f = fit(fp_high(1:ss:end),f_high(1:ss:end),'poly1');
	L{2} = strcat('m=',oval(f.p1));
	legend(l,L,'Location','southeast')
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

ub.C = 1; ub.act = 0; ub.sigma = 4;

% extra bounds
lb.act = 0; ub.act = 1;
ub.hl = 10; lb.hl = 0;
lb.n_y = 2; lb.n_z = 2;
ub.n_y = 2; ub.n_z = 2;
lb.s0 = -5; ub.s0 = 1;
ub.tau_z = 1; ub.tau_y = 1;
ub.frac = 1; lb.frac = 0;


function [] = closess(~,~)
	delete(DASim_figure)
	evalin('base','clear DASim_figure')
	evalin('base','clear DASim_plot')

end


end



