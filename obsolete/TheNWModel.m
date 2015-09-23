% The Nagel-Wilson Model
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% internal housekeeping: determine if being called by publish or not
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end

%% The Nagel-Wilson Model
% In this document, we investigate the model of Nagel and Wilson and ask: 
%
% # what the response of the model is to pulses of odor of different amplitude 
% # what the response of the model to Carlotta's random flickering stimulus is
% # if this model, for these chosen parameters, shows some sort of fast gain control

% ########   #######   ######  ######## 
% ##     ## ##     ## ##    ## ##       
% ##     ## ##     ## ##       ##       
% ##     ## ##     ##  ######  ######   
% ##     ## ##     ##       ## ##       
% ##     ## ##     ## ##    ## ##       
% ########   #######   ######  ######## 


% ########  ########  ######  ########   #######  ##    ##  ######  ######## 
% ##     ## ##       ##    ## ##     ## ##     ## ###   ## ##    ## ##       
% ##     ## ##       ##       ##     ## ##     ## ####  ## ##       ##       
% ########  ######    ######  ########  ##     ## ## ## ##  ######  ######   
% ##   ##   ##             ## ##        ##     ## ##  ####       ## ##       
% ##    ##  ##       ##    ## ##        ##     ## ##   ### ##    ## ##       
% ##     ## ########  ######  ##         #######  ##    ##  ######  ######## 


%% Response of NW Model to different pulses of Odor. 
% The N-W model consists of many parts:
%
% # a two-state receptor-ligand binding system (4 ODEs)
% # a diffusable factor that affects channel opening probabilites (2 ODEs)
% # a simple point neuron model that converts channel opening probabilities to a voltage (1 ODE)
% # a "spike" filter that converts this voltage into a firing rate (a 1-D filter)

%%
% These are the parameters used for the model:

p.ka = .001; 
p.sa = 10; % 1/s
p.kb = 1;
p.sb = 100; % 1/s
p.theta = 100; % scaling factor
p.ko = 500;
p.kc = 10;
p.a = .8; % alpha, in the diffusible factor eq. (1/s)
p.b = .6; % beta, in the diffusible factor eq. (1/s)
p.tau = .001;
p.V0 = -70;
p.R = 1;
p.Ec= 10;
p.stim_scale = 1;
p.f0 = 10; % mean output firing rate

disp(p)

%%
% We now solve the model for pulses of different heights (in practise, we vary the "stimulus scale" parameter). The following figure shows the response of the model to 1-s pulses of different heights. The model goes from a tonic response to a phasic response, whose timescale is defined by the time-scale of the "spike" filter (<20ms). The dose-response of the model looks sigmoidal.  

conc = -5:5; % in log scale
t = 0:1e-3:5;
f = NaN(11,length(t));
o = 0*t;
o(1000:2000) = 1;


for i = 1:length(conc)
	p.stim_scale = conc(i);
	[f_this,T2] = SolveNWModel(t,o,p,[]);
	f(i,:) = interp1(T2,f_this,t);
	f(i,:) = f(i,:) - min(f(i,1:100));
end

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
colors = jet(length(conc));
subplot(1,2,1), hold on
for i = 1:length(conc)
	plot(t,f(i,:),'Color',colors(i,:))
end
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')

subplot(1,2,2), hold on
temp =max(f(:,1:5000)');
plot(conc,temp,'k')
for i = 1:length(conc)
	scatter(conc(i),temp(i),532,colors(i,:),'filled')
end
xlabel('Log Concentration')
ylabel('Peak Firing rate (Hz)')

PrettyFig;


if being_published
	snapnow;
	delete(gcf)
end

return

%         ########     ###    ##    ## ########      ######  ######## #### ##     ## 
%         ##     ##   ## ##   ###   ## ##     ##    ##    ##    ##     ##  ###   ### 
%         ##     ##  ##   ##  ####  ## ##     ##    ##          ##     ##  #### #### 
%         ########  ##     ## ## ## ## ##     ##     ######     ##     ##  ## ### ## 
%         ##   ##   ######### ##  #### ##     ##          ##    ##     ##  ##     ## 
%         ##    ##  ##     ## ##   ### ##     ##    ##    ##    ##     ##  ##     ## 
%         ##     ## ##     ## ##    ## ########      ######     ##    #### ##     ## 
        
%% Response of NW Model to Carlotta's flickering stimulus
% We now look at how this model response to Carlotta's randomly flickering stimulus. To get some reasonable output from the model, we modify the parameters as follows:

p.ka = .0010;
p.f0 = 15;
p.theta = 99.1511;
p.sb = 800;


disp(p)

%%
% The only changes are a very small increase to $k_{a}$, to prevent a runaway reaction, and adjusting the arbitrary constant we add to the firing rate output. We also speed up the dynamics of the diffusible factor $D$, as we expect the adaptation from this happens on a timescale similar to that of the response dynamics. 

load('NWtestdata.mat','o','t');
[f_this,T] = SolveNWModel(t,o,p,[]);
f = interp1(T,f_this,t);
f(f<0) = 0;

% crop to fit the data
td=4;
data(td).NWORN = interp1(t,f,data(td).time);

figure('outerposition',[0 0 1000 800],'PaperUnits','points','PaperSize',[1000 800]); hold on
subplot(2,1,1), hold on
plot(data(td).time,data(td).PID,'k')
xlabel('Time (s)')
ylabel('Stimulus (V)')
set(gca,'XLim',[15 20])
subplot(2,1,2), hold on
plot(data(td).time,data(td).ORN,'k')
plot(data(td).time,data(td).NWORN,'r')
xlabel('Time (s)')
ylabel('NW ORN Firing Rate (Hz)')
set(gca,'XLim',[15 20])

PrettyFig;


if being_published
	snapnow;
	delete(gcf)
end




