% findFudgeFactor.m
% finds the best fudge factor so that when we switch variances, the mean remains the same
% 
% created by Srinivas Gorur-Shandilya at 1:25 , 09 December 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

fudge_factor_range = (0.7:0.025:.9);

% first, use manual control to make sure you have depleted the headspace
% nicely. 

alldata = struct;

for i = 1:length(fudge_factor_range)
    fudge_factor = fudge_factor_range(i);
    
    % make control paradigm with this
    ControlParadigm = makeVarianceSwitching(fudge_factor);
    
    % now run the experiment
    alldata(i).data = kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',[2 2],'w',10000);
    
end

% now make a plot of the mean and the variance in each epoch as a function of the fudge factor
mean_ratio = [];
fudge_factor = [];
for i = 1:length(alldata)
	PID = alldata(i).data(2).PID;
	for j = 1:width(PID)
		% reshape into 5-second blocks
		temp = (reshape(PID(j,:)',5e4,length(PID)/5e4));
		temp(:,1) = [];
		temp(:,end) = [];

		mean_ratio = [mean_ratio mean(temp(:,1:2:end))./mean(temp(:,2:2:end))];
		fudge_factor = [fudge_factor fudge_factor_range(i)*ones(1,width(temp)/2)];
	end
end

figure, hold on
plot(mean_ratio,fudge_factor,'+')
ff = fit(mean_ratio(:),fudge_factor(:),'poly1');
plot([min(mean_ratio) max(mean_ratio)],ff([min(mean_ratio) max(mean_ratio)]),'r')

% show what the value is at mean_ratio 1
ff(1)
