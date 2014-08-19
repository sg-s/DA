% PulsePeakAnalysis.m
% performs Pulse-peak-Analysis of gain for a flickering stimulus dataset.
% Similar to gain analysis, it analyses how the gain (defined by the ratio of the response to the stimulus) varies as a function of the stimulus in some recent history length. 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [m,fitq] = PulsePeakAnalysis(x,history_lengths,example_history_length,ph)

m = NaN*history_lengths;
fitq = NaN*history_lengths;

% set defaults
marker_size=10;
marker_size2=24;
font_size=20;
plotid=[1 2];


if isempty(ph)
	figure; hold on;
	ph(1) = gca;
	figure; hold on;
	ph(2) = gca;
end

% unpack data
f = x.response(:);
stimulus = x.stimulus(:);
t = x.time(:);
valve = x.valve(:);


% figure out sampling rate
dt = mean(diff(t));
hl = round(history_lengths/dt);


[ons,offs] = ComputeOnsOffs(valve);
ORNPeaks = NaN(1,length(ons));
PIDPeaks = NaN(1,length(ons));
Gain =  NaN(1,length(ons));
MeanStimulus = NaN(1,length(ons));

w = floor(example_history_length/dt);
s = find(ons>w,1,'first')+1;

for i = s:length(ons)
	ORNPeaks(i)=max(f(ons(i):offs(i)));
	PIDPeaks(i)=max(stimulus(ons(i):offs(i)));
	MeanStimulus(i) = mean(stimulus(ons(i)-w:ons(i)));
end

Gain = ORNPeaks./PIDPeaks;
ORNPeaks(1:s-1) = [];
PIDPeaks(1:s-1) = [];
Gain(1:s-1) = [];
MeanStimulus(1:s-1) = [];

% fit a line
[ffit, gof] = fit(MeanStimulus(:),Gain(:),'Poly1');
x = sort(unique(MeanStimulus));
y = ffit(x);


plot(ph(1),MeanStimulus,Gain,'k.','MarkerSize',marker_size2)
plot(ph(1),x,y,'r')
ylabel(ph(1),'Gain (response peak/stimulus peak)')
xlabel(ph(1),'Mean Stimulus in history window')
title(ph(1),strcat('History Length:',oval(example_history_length,2),'s'))


% now do it for all history lengths requested. 
for j = 1:length(history_lengths)
	w = hl(j);
	s = find(ons>w,1,'first')+1;

	ORNPeaks = NaN(1,length(ons));
	PIDPeaks = NaN(1,length(ons));
	Gain =  NaN(1,length(ons));
	MeanStimulus = NaN(1,length(ons));

	for i = s:length(ons)
		ORNPeaks(i)=max(f(ons(i):offs(i)));
		PIDPeaks(i)=max(stimulus(ons(i):offs(i)));
		MeanStimulus(i) = mean(stimulus(ons(i)-w:ons(i)));
	end

	Gain = ORNPeaks./PIDPeaks;
	ORNPeaks(1:s-1) = [];
	PIDPeaks(1:s-1) = [];
	Gain(1:s-1) = [];
	MeanStimulus(1:s-1) = [];


	% fit a line
	[ffit, gof] = fit(MeanStimulus(:),Gain(:),'Poly1');

	m(j) = ffit.p1;
	fitq(j) = gof.rsquare;


end

% make the second plot
hx=plotyy(ph(2),history_lengths,m,history_lengths,fitq);
ylabel(hx(1),'Slope gain/Stimulus')
ylabel(hx(2),'Quality of Fit (r-square)')
xlabel(ph(2),'History Length (s)')
