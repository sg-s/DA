% makes a phase plot of offset vs. slope for experiments where we use light and odor, one flickering and one constant. 
% 
% created by Srinivas Gorur-Shandilya at 10:34 , 01 May 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [] = GainPhasePlot(alldata,c)


for i = 1:length(alldata)
	gain = NaN(length(unique(alldata(i).paradigm)),1);
	gain_err = gain;
	offset = gain;
	offset_err = gain; 


	% first grab some data about the no background case (so either odor flicker alone or light flicker alone)
	resp0= alldata(i).resp(:,alldata(i).paradigm==1);
	stim0= alldata(i).stim(:,alldata(i).paradigm==1);
	if width(resp0) > 1
		resp0 = mean2(resp0);
		stim0 = mean2(stim0);
	end

	
	for j = 1:length(alldata(i).paradigm) % for each paradigm
		temp_gain = [];
		temp_offset = [];
		for k = find(alldata(i).paradigm==j)
			% calculate stimulus variations

			y = alldata(i).resp(:,k);
			y = y(:);
			x = resp0(:); 
			ff = fit(x,y,'poly1');
			temp_gain = [temp_gain ff.p1];
			temp_offset = [temp_offset mean(y)];


		end
		offset(j) = mean(temp_offset);
		stim_err (j) = std(temp_offset)/(sqrt(length(temp_offset)));

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

	plot(gain,offset,'-+','Color',c)


end	

