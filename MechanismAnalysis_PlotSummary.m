% 
% 
% created by Srinivas Gorur-Shandilya at 10:34 , 01 May 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [] = MechanismAnalysis_PlotSummary(alldata)

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on

offset = 3e4; % first x seconds nuked

for i = 1:length(alldata)
	subplot(1,length(alldata),i), hold on
	gain = NaN(length(alldata(i).paradigm),1);
	gain_err = gain;
	mean_stim = gain;
	stim_err = gain; 

	resp0= alldata(i).resp(offset:end,alldata(i).paradigm==1);
	stim0= alldata(i).stim(offset:end,alldata(i).paradigm==1);
	if width(resp0) > 1
		resp0 = mean2(resp0);
		stim0 = mean2(stim0);
	end

	
	for j = 1:length(alldata(i).paradigm) % for each paradigm
		temp_gain = [];
		temp_stim = [];
		for k = find(alldata(i).paradigm==j)
			% calculate stimulus variations
			temp_stim = [temp_stim mean(alldata(i).stim(offset:end,k))];

			y = alldata(i).resp(offset:end,k);
			y = y(:);
			x = resp0(:); 
			ff = fit(x,y,'poly1');
			temp_gain = [temp_gain ff.p1];

		end
		mean_stim(j) = mean(temp_stim);
		stim_err (j) = std(temp_stim)/(sqrt(length(temp_stim)));

		gain(j) = mean(temp_gain);
		gain_err(j) = std(temp_gain)/(sqrt(length(temp_gain)));

		

	end

	xt = {};
	for j = 1:length(alldata(i).ParadigmNames)
		t = alldata(i).ParadigmNames{j};
		a = strfind(t,'+');
		z = strfind(t,'V');
		xt{j} = t(a+1:z-1);
	end

	mean_stim = mean_stim/mean(stim0);
	stim_err = stim_err/mean(stim0);
	errorbar(mean_stim,stim_err,'b')
	errorbar(gain,gain_err,'r')

	ylabel('Fold Change')
	xlabel('Supplemental Light (V)')
	legend('Mean Stimulus','Gain')
	set(gca,'XTick',[1:length(alldata(i).ParadigmNames)],'XTickLabel',xt,'XLim',[0 6],'YLim',[.2 1.2])

	

end	

