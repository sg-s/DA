% InteractiveGainAnalysis.m
% this is an interactive version of the gain analysis function, uses Manipulate.m
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [] = InteractiveGainAnalysis(p)

% get the data to work with from the base workspace
data = evalin('base', 'IGA_data');

% check if the figure has been created
IGA_figure = [];
try IGA_figure = evalin('base', 'IGA_figure');
	IGA_plot = evalin('base', 'IGA_plot');
catch
	IGA_figure = figure('Units','pixels','Position',[10 100 1149 700],'Name','Interactive Gain Analysis','IntegerHandle','off','CloseRequestFcn',@closess,'Color','w'); hold on
	clf(IGA_figure);
	IGA_plot(1)=axes('Units','pixels','Position',[52.5 490 612.5 175]);
	plot(data.time,data.stimulus,'k')
	IGA_plot(2)=axes('Units','pixels','Position',[52.5 280 612.5 175]);
	plot(data.time,data.response,'k')
	IGA_plot(3)=axes('Units','pixels','Position',[52.5 52.5 612.5 192.5]);
	plot(data.time,data.prediction,'k')
	linkaxes(IGA_plot(1:3),'x');
	linkaxes(IGA_plot(2:3),'y');
	set(IGA_plot(1),'box','off','XLim',[min(data.time) max(data.time)])
	set(IGA_plot(2),'box','off','XLim',[min(data.time) max(data.time)])
	set(IGA_plot(3),'box','off','XLim',[min(data.time) max(data.time)])
	IGA_plot(4)=axes('Units','pixels','Position',[717.5 332.5 385 332.5]);
	IGA_plot(5)=axes('Units','pixels','Position',[717.5 52.5 385 210]);
	ylabel(IGA_plot(1),'Stimulus')
	ylabel(IGA_plot(2),'Response')
	ylabel(IGA_plot(3),'Prediction')

	% save to the base workspace
	assignin('base','IGA_figure',IGA_figure)
	assignin('base','IGA_plot',IGA_plot)

	% clean up data
	rm_this = unique([find(isnan(data.response)) find(isnan(data.prediction)) find(isnan(data.stimulus))]);
	data.response(rm_this) = [];
	data.prediction(rm_this) = [];
	data.stimulus(rm_this) = [];
	data.time(rm_this) = [];
end



% say what the parameters are for getModelParameters
p.frac;
[~,~,~,hl_min]=FindCorrelationTime(data.stimulus);
dt = mean(diff(data.time));
hl_min = (hl_min*dt);
p.hl_max; % in s
history_lengths = logspace(log10(hl_min),log10(p.hl_max),30);
p.ehl;
lb.ehl = 1;
ub.ehl = 30;
ub.frac = 1;
lb.frac = 0;
ub.hl_max = 10;
lb.hl_max = 1; 

% pick the example history length
example_history_length = history_lengths(round(p.ehl));

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

	cla(IGA_plot(4))
	cla(IGA_plot(5))
	% do the gain analysis. 
	ph(3:4) = IGA_plot(4:5);
	p_Linear = NaN(2,30);
	
	data.frac = p.frac;
	GainAnalysis4(data,history_lengths,example_history_length,ph,p_Linear);
	ylabel(IGA_plot(5),'Gain')
	set(IGA_plot(5),'XScale','log')

	% Indicate regions of high and low stimulus on the stimulus
	hl = floor(example_history_length/dt);
	shat = ComputeSmoothedStimulus(data.stimulus,hl);

	n = floor(p.frac*length(data.response));
	shat(1:hl) = Inf; % the initial segment where we can't estimate shat is excluded
	shat(isnan(shat)) = Inf;
	shat(isnan(data.response)) = Inf;
	[~, t_low] = sort(shat,'ascend');
	t_low = t_low(1:n); % this is an index
	t_low = data.time(t_low); % t_low is now a time. 
	 
	shat = ComputeSmoothedStimulus(data.stimulus,hl);
	shat(1:hl) = -Inf;
	shat(isinf(shat)) = -Inf;
	shat(isnan(data.response)) = -Inf;
	[~, t_high] = sort(shat,'descend');
	t_high = t_high(1:n);
	t_high  = data.time(t_high);

	ypos = max(data.stimulus) + 0.1*std(data.stimulus);
	cla(IGA_plot(1))
	hold(IGA_plot(1),'on')
	plot(IGA_plot(1),data.time,data.stimulus,'k')
	plot(IGA_plot(1),t_low,ypos+0*t_low,'g.')
	plot(IGA_plot(1),t_high,ypos+0*t_low,'r.')


end

function [] = closess(~,~)
	delete(IGA_figure)
	evalin('base','clear IGA_figure')
	evalin('base','clear IGA_plot')

end

end