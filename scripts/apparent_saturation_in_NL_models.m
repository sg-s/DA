
%% Apparent saturation in NL models
% In this document, I look at a property of NL models that allow them to code for the duration of brief pulses of signal in the amplitude of the response. This effect is pronounced when signals are close to saturating, and when the duration of the signals is less than the timescale of the filter in the NL model. 

pHeader;

% these parameters come from the model fit to MSG data

k_D = 0.0968;

clear p
p.tau1 = 30;
p.tau2 = 100;
p.     A = 0.7000;
p.     n = 2;

% generate the filter
K = filter_gamma2(1:800,p);

pulse_durations = round(logspace(0,3,100));
pulse_amplitudes = logspace(-3,3,100); % in units of k_D

R = NaN(length(pulse_durations),length(pulse_amplitudes));

for i = 1:length(pulse_durations)
	for j = 1:length(pulse_amplitudes)
		% make the stimulus
		S = zeros(5e3,1);
		S(1e3:1e3+pulse_durations(i)) = pulse_amplitudes(j)*k_D;

		% pass through nonlinearity 
		a = 1./(1+k_D./S);

		% filter it 
		R(i,j) = max(filter(K,1,a));
	end
end

%%
% in the following figure, I generate responses from a NL model and quantify the response of the model as the peak response for a pulse with a given duration and amplitude. (a) Peak response of NL model as a function of pulse amplitude, for various pulse durations. Pulse durations are in units of timescale of the filter $\tau_K$, and pulse amplitudes are in units of $k_D$. Note that it looks like the NL model has an apparent saturation below its true saturation for short pulses. (b) Peak response as a function of pulse duration, for pulses of varying amplitudes. Note that the amplitude of the NL model response scales with the duration of the pulse. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on

% R vs amplitude, colours are various durations 
subplot(1,2,1); hold on
c = parula(length(pulse_durations));
clear l L
j = 1;
for i = 1:10:length(pulse_durations)
	l(j) = plot(pulse_amplitudes,R(i,:),'Color',c(i,:));
	L{j} = ['T/\tau_{K} = ' oval(pulse_durations(i)/p.tau2)];
	j = j + 1;
end
legend(l,L,'Location','northwest')
set(gca,'XScale','log','XLim',[min(pulse_amplitudes) max(pulse_amplitudes)],'XTick',[1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3])
xlabel('S/k_D')
ylabel('R')

% R vs duration, colours are various amplitudes 
subplot(1,2,2); hold on
c = jet(length(pulse_amplitudes));
clear l L
j = 1;
for i = 35:5:65
	l(j) = plot(pulse_durations/p.tau2,R(:,i),'Color',c(i,:));
	L{j} = ['S/K_{D} = ' oval(pulse_amplitudes(i))];
	j = j + 1;
end
legend(l,L,'Location','southeast')
set(gca,'XScale','log','XLim',[.01 10],'XTick',[.01 .1 1 10],'YScale','log')
xlabel('T/\tau_{K}')
ylabel('R')

prettyFig()

labelFigure

if being_published	
	snapnow	
	delete(gcf)
end

%%
% Can this phenomenon be used in the actual neuron? Well, there is a complication: if both the amplitude and duration of the stimulus are broadly distributed, and can vary quickly, a given response can correspond to a short signal with a large amplitude, or a long signal with a small amplitude. There seems to be no way to distinguish these two cases. So this coding scheme would only be useful in scenarios where either the pulse amplitude or the pulse duration is fixed. 

%%
% If you have extremely sensitive receptors, then you are in the regime where the effective pulse amplitude is fixed, as every pulse saturates your receptors. 

%% Version Info
%
pFooter;


