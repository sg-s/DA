%% Does Spike Frequency Adaptation imply Dynamical Adaptation? 
% Liu and Wang's (papers2://publication/uuid/54F8985C-DC25-42D9-AD8B-805EBCD7475F) model of spike frequency adaptation shows some qualitative similarity to ORN response dynamics. Can the mechanics of spike frequency adaptation also give rise to a phenomenon like fast dynamical adaptation? 

%% Response properties of Model Neuron 
% In this plot, a Liu-Wang model neuron is simulated using parameters very close to what they use in the paper. The parameters are:

p.Cm=0.5000;
p.gL=0.0250;
p.tau_Ca; 1;
p.C=0.3800;
p.gAHP=1.3500;
p.Vreset= -60;
p.Vth= -55;
p.Vrest= -70;
p.A= 666;
p.Vk= -80;

disp(p)

%%
% The figure below shows the response of this model neuron to increasing pulses of stimulus. 
time = 1e-4:1e-4:5;
stim = zeros(length(time),10);
f = zeros(length(time),10);
for i = 1:10
	clear temp
	stim(10000:30000,i) = i*0.1;
	nrep = 50;
	st = zeros(length(time),nrep);
	for j = 1:nrep
		[~, ~,st(:,j)] = XJWNeuronEuler(time,stim(:,i)+0.3*randn(length(time),1),p);
		temp(:,j)=spiketimes2f(st(:,j),time,0.003,0.03);
	end
	temp = mean(temp');
	t=0:0.003:5;
	f(:,i) = interp1(t,temp,time);

end



figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
plot(time,stim)
xlabel('Time (s)')
ylabel('Stimulus (a.u.)')
PrettyFig;

subplot(1,2,2), hold on
plot(time,f)
xlabel('Time (s)')
ylabel('ORN Response (Hz)')
PrettyFig;

%% Response to exponentiated Gaussian stimuli
% Now we simulate some flickering stimuli by exponentiated Gaussian noise and use that to feed the neuron. The following figure shows the output of the model neuron, together with a linear fit to the model output. 

time = 1e-4:1e-4:30;
stim = exp(randn(length(time),1));
stim = filter(ones(100,1)/100,1,stim);
nrep = 10;
f = zeros(length(time),nrep);
st = zeros(length(time),nrep);
t=0.003:0.003:30;
temp = zeros(length(t),nrep);
for j = 1:nrep
	[~, ~,st(:,j)] = XJWNeuronEuler(time,stim+0.1*randn(length(time),1),p);
	temp(:,j)=spiketimes2f(st(:,j),time,0.003,0.03);
end
f = mean(temp');
stim = interp1(time,stim,t);

f(t<5) = [];
stim(t<5) = [];
t(t<5)=[];


% fit a linear model
[K,~,filtertime] = FindBestFilter(stim,f);
fp = mean(f)+convolve(t,stim,K,filtertime);

figure('outerposition',[0 0 600 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(t,f,'k'), hold on
plot(t,fp,'r')
legend Model LinearPrediction
xlabel('Time (s)')
ylabel('ORN Response (Hz)')
PrettyFig;
set(gca,'XLim',[20 25])

%% Does this show the signature of fast adaptation?
% Does the (feedback) mechanism of spike frequency adaptation allow the neuron to modulate its gain on a fast time scale, like in real ORNs? We perform a linear gain analysis as before. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on

s = 300; % when we start for the gain analysis
z = length(f); % where we end
example_history_length = 1.02;
history_lengths = [0 0.003 0.03 0.3 0.6 0.711 1.02 1.2 1.5 2.1];

clear x
x.response = f(s:z);
x.prediction = fp(s:z);
x.stimulus = stim(s:z);
x.time = t(s:z);
x.filter_length = 201;

redo_bootstrap = 0;
if redo_bootstrap
	data(td).LinearFit_p = GainAnalysis3(x,history_lengths,example_history_length,ph);
else
	GainAnalysis3(x,history_lengths,example_history_length,ph,Inf*history_lengths);
end
clear x
hold off
set(ph(4),'YLim',[0.7 1.4],'XLim',[-0.01 max(history_lengths)])