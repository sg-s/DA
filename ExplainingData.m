% ExplainingData
% this m file attempts to fit models to all the data we present, and explain the observed phenemnology 
% 
% created by Srinivas Gorur-Shandilya at 12:52 , 18 March 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

%              ###     ######   ######  ######## ##     ## ########  ##       ######## 
%             ## ##   ##    ## ##    ## ##       ###   ### ##     ## ##       ##       
%            ##   ##  ##       ##       ##       #### #### ##     ## ##       ##       
%           ##     ##  ######   ######  ######   ## ### ## ########  ##       ######   
%           #########       ##       ## ##       ##     ## ##     ## ##       ##       
%           ##     ## ##    ## ##    ## ##       ##     ## ##     ## ##       ##       
%           ##     ##  ######   ######  ######## ##     ## ########  ######## ######## 

%           ########     ###    ########    ###    
%           ##     ##   ## ##      ##      ## ##   
%           ##     ##  ##   ##     ##     ##   ##  
%           ##     ## ##     ##    ##    ##     ## 
%           ##     ## #########    ##    ######### 
%           ##     ## ##     ##    ##    ##     ## 
%           ########  ##     ##    ##    ##     ## 


load('MeanShiftedGaussians.mat')

% shorten paradigm names by throwing out 'MFC'
short_paradigm_names = paradigm_names;
for i = 1:length(paradigm_names)
	short_paradigm_names{i} = paradigm_names{i}(strfind(paradigm_names{i},'-')+1:end);
end

detrended_data = cache('detrended_data'); 


%       ########     ###               ##       ##     ##  ######   ######   
%       ##     ##   ## ##               ##      ###   ### ##    ## ##    ##  
%       ##     ##  ##   ##               ##     #### #### ##       ##        
%       ##     ## ##     ##    #######    ##    ## ### ##  ######  ##   #### 
%       ##     ## #########              ##     ##     ##       ## ##    ##  
%       ##     ## ##     ##             ##      ##     ## ##    ## ##    ##  
%       ########  ##     ##            ##       ##     ##  ######   ######   

% Fit DA Model to Mean Shifted Gaussians 

load('DA_Fit_to_MSG.mat')

for i = 1:3
	% prep data
	clear d
	c = 1;
	for j = 1:8
		if ~isempty(detrended_data(j,i).stim)
			d(c).stimulus = detrended_data(j,i).stim;
			d(c).response = detrended_data(j,i).resp;
			d(c).response(1:1e3) = NaN;
			c=c+1;
		end
	end

	% optimise
	for k = 1:5
		p(i) = FitModel2Data(@DAModelv2,d,p(i));
	end

	save('DA_Fit_to_MSG.mat','p','-append')
end


%%
% The following plot shows the gain of the ORN vs stimulus, with the best-fit DA model overlaid. Data is in black, and the DA model fits are in red. The data is segmented into experimental "days", so that the same set of neurons appears in all stimulus means. The DA model is also fit in this fashion, so we get three sets of DA model fits. 

figure('outerposition',[0 0 600 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
clear l
% plot the gain of the data
for i = 1:3
	stim = [];
	resp = [];
	for j = 1:8
		if ~isempty(detrended_data(j,i).stim)
			stim = [stim; detrended_data(j,i).stim];
			resp = [resp detrended_data(j,i).resp];
		end
	end
	g = EstimateGain2(stim,resp);
	l(1)=plot(mean(stim'),g,'k+');
end

% plot the gain of the DA model fits
for i = 1:3
	stim = [];
	resp = [];
	for j = 1:8
		if ~isempty(detrended_data(j,i).stim)
			stim = [stim; detrended_data(j,i).stim];
			this_resp = DAModelv2(detrended_data(j,i).stim,p(i));
			this_resp(1:1e3) = NaN;
			resp = [resp; this_resp];
		end
	end
	g = EstimateGain2(stim,resp);
	l(2)=plot(mean(stim'),g,'r+');
end

set(gca,'XScale','log','YScale','log')
xlabel('Mean Stimulus (V)')
ylabel('Gain (Hz/V)')
legend(l,{'Data','DA Model'})

PrettyFig('plw=1.5;','lw=1.5;','fs=14;')

if being_published
	snapnow
	delete(gcf)
end



