
%% The Nagel-Wilson Model
% In this document we try to make sense of the Nagel Wilson Model and see if it explains gain scaling and LFP slowdown.

pHeader;

% core parameters 
T = 20e3; % total length
pulse_on = 15e3;
pulse_off = 16e3;

%% 1: Using the Nagel-Wilson Model as is
% In this section, we use the Nagel-Wilson model as described in the paper. This has 6 ODEs, and the output of the model is the channel open probability. In the following figure, I plot the response of this model to a test pulse of stimulus. 

clear p
p.theta = 100;
p.   ka = 0.001;
p.   kb = 1;
p.   sa = 10;
p.   sb = 100;
p.   k0 = 500;
p.    A = 0.8000;
p.    B = 0.6000;
p.   ko = 500;
p.   kc = 10;

S = .1+zeros(T,1);
S(pulse_on:pulse_off) = 1;
[C,D,R] = NagelWilsonIntegrate(S,p);

figure('outerposition',[0 0 1200 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
subplot(2,2,1), hold on
plot(S)
ylabel('Stimulus')
subplot(2,2,2), hold on
plot(R)
legend({'R','R*','OR','OR*'},'Location','northwest')
ylabel('Receptor Species')
subplot(2,2,3), hold on
plot(C)
ylabel('C_{open}')
subplot(2,2,4), hold on
plot(D)
ylabel('Diffusible factor')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% We see that there is a fatal flaw in the model: there is no conservation equation for the four reaction species. Thus, they are allowed to grow without bound, which is obviously wrong. 

%% 2. Nagel Wilson Model with conservation equation
% In this section, we add a conservation equation to the model:
% 
% $$ OR=R_{total}-R-R*-OR*$$
%

%%
% Doing so eliminates one of the ODEs. The ODEs now are:
%
% $$ \dot{R}=-R(k_{a}s_{a}+Ok_{b}s_{b})+R*s_{a}+(R_{total}-R-R*-OR*)s_{b} $$
%

%%
% and so on.

%%
% In the following figure, we plot the responses to a pulse of odorant.

S = .1+zeros(T,1);
S(pulse_on:pulse_off) = .2;
[C,D,R] = NagelWilsonModelReduced(S,p);

figure('outerposition',[0 0 1200 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
subplot(2,2,1), hold on
plot(S)
set(gca,'XLim',[length(S)-10e3 length(S)])
ylabel('Stimulus')
subplot(2,2,2), hold on
plot(R)
set(gca,'XLim',[length(S)-10e3 length(S)])
legend({'R','R*','OR','OR*'})
ylabel('Receptor Species')
subplot(2,2,3), hold on
plot(C)
set(gca,'XLim',[length(S)-10e3 length(S)])
ylabel('C_{open}')
subplot(2,2,4), hold on
plot(D)
set(gca,'XLim',[length(S)-10e3 length(S)])
ylabel('Diffusible factor')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% OK, this is somewhat reasonable. Now, we check to see if this model predicts a slowdown in response with increased stimulus. We do this by stimulating the model with pulses on top of increasing backgrounds. We also check to see if the responses follow Weber's Law.

clear ax
figure('outerposition',[0 0 1200 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
for i = 1:6
  ax(i) = subplot(2,3,i); hold on
end

background_levels = logspace(-2,0,10);
c = parula(length(background_levels)+1);

receptor_gain = NaN*background_levels;
channel_gain = NaN*background_levels;
for i = 1:length(background_levels)
	S = background_levels(i) + zeros(T,1);
	S(pulse_on:pulse_off) = 2*background_levels(i);
	plot(ax(1),S,'Color',c(i,:));

	[C,D,R] = NagelWilsonModelReduced(S,p);
	plot(ax(2),R(:,2)+R(:,4),'Color',c(i,:))

	plot(ax(3),C,'Color',c(i,:))

	plot(ax(4),D,'Color',c(i,:))

	% compute gain at receptors 
	y = R(:,2) + R(:,4);
	delta_r = max(y(pulse_on:pulse_off)) - y(pulse_on-1);
	g = delta_r/(max(S) - min(S));
	receptor_gain(i) = g;
	plot(ax(5),min(S),g,'+','Color',c(i,:));

	% plot total gain
	delta_r = max(C(pulse_on:pulse_off)) - C(pulse_on-1);
	g = delta_r/(max(S) - min(S));
	channel_gain(i) = g;
	plot(ax(6),min(S),g,'+','Color',c(i,:));
end

for i = 1:4
	set(ax(i),'XLim',[length(S) - 10e3 length(S)])
end

% draw a weber's fit to receptor gain
x = background_levels;
y = receptor_gain;
options = fitoptions(fittype('poly1'));
options.Lower = [-1 -Inf];
options.Upper = [-1 Inf];
cf_temp = fit(log(x(:)),log(y(:)),'poly1',options);
cf = fit(x(:),y(:),'power1');
warning off
cf.a = exp(cf_temp.p2); cf.b = -1;
warning on
plot(ax(5),sort(x),cf(sort(x)),'r')

% and also for channel gain
y = channel_gain;
options = fitoptions(fittype('poly1'));
options.Lower = [-1 -Inf];
options.Upper = [-1 Inf];
cf_temp = fit(log(x(:)),log(y(:)),'poly1',options);
cf = fit(x(:),y(:),'power1');
warning off
cf.a = exp(cf_temp.p2); cf.b = -1;
warning on
plot(ax(6),sort(x),cf(sort(x)),'r')

set(ax(1),'YScale','log')
ylabel(ax(1),'Stimulus')
ylabel(ax(2),'R* + OR*')
ylabel(ax(3),'C_{open}')
ylabel(ax(4),'Diffusible factor')

title(ax(5),'Receptor Gain')
ylabel(ax(5),'\DeltaR/\DeltaS')
set(ax(5),'YScale','log','XScale','log')
xlabel(ax(5),'Background Stim')

title(ax(6),'Channel Gain')
ylabel(ax(6),'\DeltaR/\DeltaS')
set(ax(6),'YScale','log','XScale','log')
xlabel(ax(6),'Background Stim')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% This model doesn't exactly follow Weber's Law, at least for these parameters, but it's close. We note that the channels alone seem to be obeying Weber's Law. We consider only the channel opening and the diffusible factor ODEs, and assume that receptor binding is so fast that the activated receptor species are exactly some scaled version of the stimulus. Then, the system collapses to just the last two ODEs in their model. This simplification leaves just four parameters, and we can show that this reduced model can show Weber-Fechner gain scaling. 

% generate model responses
clear p
p.    r_b = 20;
p.    r_d = 10;
p.theta_b = .05;
p.theta_d = 1e-1;


clear ax
figure('outerposition',[0 0 800 700],'PaperUnits','points','PaperSize',[800 700]); hold on
for i = 1:4
	ax(i) = subplot(2,2,i); hold on
end

background_levels = logspace(-2,0,10);
c = parula(length(background_levels)+1);

all_gain = NaN*background_levels;
for i = 1:length(background_levels)
	S = background_levels(i) + zeros(T,1);
	S(pulse_on:pulse_off) = 2*background_levels(i);
	plot(ax(1),S,'Color',c(i,:));

	[R,D] = simpleReceptorModelX(S,p);
	plot(ax(2),R,'Color',c(i,:))

	plot(ax(3),D,'Color',c(i,:))

	% plot total gain
	delta_r = max(R(pulse_on:pulse_off)) - R(pulse_on-1);
	g = delta_r/(max(S) - min(S));
	all_gain(i) = g;
	plot(ax(4),min(S),g,'+','Color',c(i,:));
end

for i = 1:3
	set(ax(i),'XLim',[length(S) - 10e3 length(S)])
end

% draw a weber's fit to gain
x = background_levels;
y = all_gain;
options = fitoptions(fittype('poly1'));
options.Lower = [-1 -Inf];
options.Upper = [-1 Inf];
cf_temp = fit(log(x(:)),log(y(:)),'poly1',options);
cf = fit(x(:),y(:),'power1');
warning off
cf.a = exp(cf_temp.p2); cf.b = -1;
warning on
plot(ax(4),sort(x),cf(sort(x)),'r')


set(ax(1),'YScale','log')
ylabel(ax(1),'Stimulus')
ylabel(ax(2),'C_{open}')
ylabel(ax(3),'Diffusible factor')

title(ax(4),'Gain')
ylabel(ax(4),'\DeltaR/\DeltaS')
set(ax(4),'YScale','log','XScale','log')
xlabel(ax(4),'Background Stim')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% The model with an error
% A previous implementation of the model had an error in a sign of first term in the second ODE. Bizarrely, this model had some nice properties:

clear ax
figure('outerposition',[0 0 1200 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
for i = 1:6
  ax(i) = subplot(2,3,i); hold on
end

background_levels = logspace(-2,0,10);
c = parula(length(background_levels)+1);

receptor_gain = NaN*background_levels;
channel_gain = NaN*background_levels;
for i = 1:length(background_levels)
	S = background_levels(i) + zeros(T,1);
	S(pulse_on:pulse_off) = 2*background_levels(i);
	plot(ax(1),S,'Color',c(i,:));

	[C,D,R] = NagelWilsonModelReduced_error(S,p);
	plot(ax(2),R(:,2)+R(:,4),'Color',c(i,:))

	plot(ax(3),C,'Color',c(i,:))

	plot(ax(4),D,'Color',c(i,:))

	% compute gain at receptors 
	y = R(:,2) + R(:,4);
	delta_r = max(y(pulse_on:pulse_off)) - y(pulse_on-1);
	g = delta_r/(max(S) - min(S));
	receptor_gain(i) = g;
	plot(ax(5),min(S),g,'+','Color',c(i,:));

	% plot total gain
	delta_r = max(C(pulse_on:pulse_off)) - C(pulse_on-1);
	g = delta_r/(max(S) - min(S));
	channel_gain(i) = g;
	plot(ax(6),min(S),g,'+','Color',c(i,:));
end

for i = 1:4
	set(ax(i),'XLim',[length(S) - 10e3 length(S)])
end

% draw a weber's fit to receptor gain
x = background_levels;
y = receptor_gain;
options = fitoptions(fittype('poly1'));
options.Lower = [-1 -Inf];
options.Upper = [-1 Inf];
cf_temp = fit(log(x(:)),log(y(:)),'poly1',options);
cf = fit(x(:),y(:),'power1');
warning off
cf.a = exp(cf_temp.p2); cf.b = -1;
warning on
plot(ax(5),sort(x),cf(sort(x)),'r')

% and also for channel gain
y = channel_gain;
options = fitoptions(fittype('poly1'));
options.Lower = [-1 -Inf];
options.Upper = [-1 Inf];
cf_temp = fit(log(x(:)),log(y(:)),'poly1',options);
cf = fit(x(:),y(:),'power1');
warning off
cf.a = exp(cf_temp.p2); cf.b = -1;
warning on
plot(ax(6),sort(x),cf(sort(x)),'r')

set(ax(1),'YScale','log')
ylabel(ax(1),'Stimulus')
ylabel(ax(2),'R* + OR*')
ylabel(ax(3),'C_{open}')
ylabel(ax(4),'Diffusible factor')

title(ax(5),'Receptor Gain')
ylabel(ax(5),'\DeltaR/\DeltaS')
set(ax(5),'YScale','log','XScale','log')
xlabel(ax(5),'Background Stim')

title(ax(6),'Channel Gain')
ylabel(ax(6),'\DeltaR/\DeltaS')
set(ax(6),'YScale','log','XScale','log')
xlabel(ax(6),'Background Stim')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%% Version Info
%
pFooter;
