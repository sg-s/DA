% 
% 
% created by Srinivas Gorur-Shandilya at 1:59 , 07 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


pHeader;

% load the data
if ~exist('od','var')
	p = '/local-data/DA-paper/natural-flickering/with-lfp/ab3/';
	od = raw2ORNData(p);
	% ignore all the LFP for now
	for i = 1:length(od)
		od(i).LFP = [];
	end

	% specify interesting range
	uts = zeros(length(od(1).stimulus),1);
	uts(10e3:end-5e3) = true;
	for i = 1:length(od)
		od(i).use_this_segment = uts;
	end

	% back out all filters
	od = backOutFilters(od);

	% fit NL
	od2 = fitNL(od);
end

figure('outerposition',[0 0 1500 800],'PaperUnits','points','PaperSize',[1500 800]); hold on
clear ax
ax(1) = subplot(2,2,1); hold on
ax(2) = subplot(2,2,2); hold on
for i = 3:6
	ax(i) = subplot(2,4,2+i); hold on
end

t = (1:length(od(1).stimulus))*1e-3;
plot(ax(1),t,nanmean([od.stimulus],2),'k')
xlabel(ax(1),'Time (s)')
ylabel(ax(1),'Stimulus (V)')
set(ax(1),'XLim',[0 70])

for i = 1:length(od)
	plot(ax(2),t,nanmean([od(i).firing_rate],2),'Color',[0 0 0 0.1])
end
plot(ax(2),t,nanmean([od.firing_rate],2),'Color','k')
xlabel(ax(2),'Time (s)')
ylabel(ax(2),'Firing Rate (Hz)')
set(ax(2),'XLim',[0 70])

% show response vs. linear projection; colour by mean stimulus in recent history window 
[~,excursions] = plotExcursions(od(2),ax(3),'data','firing_rate');

% fit a NL just to the excursions
temp = ORNData;
temp.stimulus = nanmean(od(2).stimulus,2);
temp.firing_rate = nanmean(od(2).firing_rate,2);
temp.firing_projected = nanmean(od(2).firing_projected,2);
temp.use_this_segment = false(length(temp.stimulus),1);
for i = 1:length(excursions.ons)
	temp.use_this_segment(excursions.ons(i):excursions.offs(i)) = true;
end
[temp2,ff] = fitNL(temp);

% show this best-fit NL on this plot
x = temp.firing_projected(temp.use_this_segment,:); x2 = x;
x2 = x2 - min(x2); x2 = x2/max(x2);
plot(ax(3),sort(x),max(temp.firing_rate(temp.use_this_segment,:))*ff(sort(x2)),'r')
xlabel(ax(3),'Proj. Stimulus (V)')
l = plot(ax(3),NaN,NaN);
linear_filter_r2 = rsquare(temp.firing_projected(temp.use_this_segment),temp.firing_rate(temp.use_this_segment));
legend(l,['r^2 = ' oval(linear_filter_r2)],'Location','southeast')

% show response vs. LN model; colour by mean stimulus in recent history window 
uts = false(length(temp.stimulus),1);
uts(10e3:end-5e3) = true;
temp2.use_this_segment = uts;
[~,excursions] = plotExcursions(temp2,ax(4),'data','firing_rate');
xlabel(ax(4),'LN Model Prediction (Hz)')
l = plot(ax(4),NaN,NaN);
LN_model_r2 = rsquare(temp2.firing_projected(temp.use_this_segment),temp2.firing_rate(temp.use_this_segment));
legend(l,['r^2 = ' oval(LN_model_r2)],'Location','southeast')

% fit a DA Model to the excursions
clear p
p.   s0 = -0.1164;
p.  n_z = 10.6250;
p.tau_z = 19.7499;
p.  n_y = 10.6250;
p.tau_y = 4.6377;
p.    C = 0.5848;
p.    A = 709.4439;
p.    B = 12.0094;
S = temp2.stimulus;
[R,y,z,Ky,Kz] = DAModelv2(S,p);
temp3 = temp2;
temp3.firing_projected = R;
plotExcursions(temp3,ax(5),'data','firing_rate');
l = plot(ax(5),NaN,NaN);
DA_model_r2 = rsquare(R(temp.use_this_segment),temp.firing_rate(temp.use_this_segment));
legend(l,['r^2 = ' oval(DA_model_r2)],'Location','southeast');
xlabel(ax(5),'DA Model Prediction (Hz)')

% now show that this is the key timescale of gain control
tau_gain = round(logspace(log10(50),4,50));
r2 = NaN*tau_gain;
for i = 1:length(tau_gain)
	p.tau_z = tau_gain(i)/p.n_z;
	R = DAModelv2(S,p);
	r2(i) = rsquare(R(temp.use_this_segment),temp2.firing_rate(temp.use_this_segment));
end

% convert into fraction remaining variance explained
r2 = (r2 - linear_filter_r2)./(1-linear_filter_r2);


plot(ax(6),tau_gain,r2,'k+')
set(ax(6),'XScale','log','YLim',[0 1])
xlabel(ax(6),'Timescale of gain control (ms)')
ylabel(ax(6),['Remaining variance' char(10) 'explained by DA model'])

% show where the LN model is on this plot
plot(ax(6),tau_gain,(LN_model_r2-linear_filter_r2)/(1-linear_filter_r2)*(1+0*tau_gain),'r')


prettyFig('fs=18;','FixLogX=true;')

if being_published	
	snapnow	
	delete(gcf)
end



%% Version Info
%
pFooter;


