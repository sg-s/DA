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

load('MSG_per_neuron.mat','MSG_data')


%       ########     ###               ##       ##     ##  ######   ######   
%       ##     ##   ## ##               ##      ###   ### ##    ## ##    ##  
%       ##     ##  ##   ##               ##     #### #### ##       ##        
%       ##     ## ##     ##    #######    ##    ## ### ##  ######  ##   #### 
%       ##     ## #########              ##     ##     ##       ## ##    ##  
%       ##     ## ##     ##             ##      ##     ## ##    ## ##    ##  
%       ########  ##     ##            ##       ##     ##  ######   ######   

% Fit DA Model to Mean Shifted Gaussians 

load('DA_Fit_to_MSG.mat')

for j = 1:13 % these many neurons
	clear d
	c = 1;
	for i = 1:8 % these many paradigms		
		if ~isempty(MSG_data(i,j).stim)
			if width(MSG_data(i,j).stim) > 1
				% we ignore neurons with only one trial
				d(c).stimulus = mean2(MSG_data(i,j).stim);
				d(c).response = mean2(MSG_data(i,j).resp);
				d(c).response(1:1e3) = NaN;
				c=c+1;
			end
		end
	end
	if exist('d','var')
		% optimise
		for k = 1:5
			p(j) = FitModel2Data(@DAModelv2,d,p(j));
		end
	end
	save('DA_Fit_to_MSG.mat','p','-append')
end





%    ##      ## ######## ########  ######## ########     ########     ###    
%    ##  ##  ## ##       ##     ## ##       ##     ##    ##     ##   ## ##   
%    ##  ##  ## ##       ##     ## ##       ##     ##    ##     ##  ##   ##  
%    ##  ##  ## ######   ########  ######   ########     ##     ## ##     ## 
%    ##  ##  ## ##       ##     ## ##       ##   ##      ##     ## ######### 
%    ##  ##  ## ##       ##     ## ##       ##    ##     ##     ## ##     ## 
%     ###  ###  ######## ########  ######## ##     ##    ########  ##     ## 



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


%   ######  ########  ######## ######## ########  ##     ## ########     ########     ###    
%  ##    ## ##     ## ##       ##       ##     ## ##     ## ##     ##    ##     ##   ## ##   
%  ##       ##     ## ##       ##       ##     ## ##     ## ##     ##    ##     ##  ##   ##  
%   ######  ########  ######   ######   ##     ## ##     ## ########     ##     ## ##     ## 
%        ## ##        ##       ##       ##     ## ##     ## ##           ##     ## ######### 
%  ##    ## ##        ##       ##       ##     ## ##     ## ##           ##     ## ##     ## 
%   ######  ##        ######## ######## ########   #######  ##           ########  ##     ## 

%%
% In this section we check if the DA model can account for observed speedup of responses with increasing stimulus mean. 

% xcorr for mean shifted gaussians
peak_loc_xcorr = NaN(length(detrended_data),3);
mean_stim = NaN(length(detrended_data),3);
clear l 
l = zeros(8,1);
for k = 1:3
	for i = 1:length(detrended_data)
		mean_stim(i,k) = mean(detrended_data(i,k).stim);
		a = detrended_data(i,k).resp - mean(detrended_data(i,k).resp);
		if ~isempty(a)
			a = a/std(a);
			b = detrended_data(i,k).stim - mean(detrended_data(i,k).stim);
			b = b/std(b);
			x = xcorr(a,b); % positive peak means a lags b
			t = 1e-3*(1:length(x));
			t = t-mean(t);
			x = x/max(x);
			[~,loc] = max(x);
			peak_loc_xcorr(i,k) = t(loc);
		end
	end
end 

% now compute xcorr for the DA model
peak_loc_xcorr_DA = NaN(length(detrended_data),3);
for i = 1:3
	for j = 1:8
		if ~isempty(detrended_data(j,i).stim)
			this_resp = DAModelv2(detrended_data(j,i).stim,p(i));
			this_resp(1:1e3) = NaN;
			detrended_data(j,i).DAFit = this_resp;
			this_stim = detrended_data(j,i).stim;
			this_stim(1:1e3) = [];
			this_resp(1:1e3) = [];
			a = this_resp - mean(this_resp);
			b = this_stim - mean(this_stim);
			a = a/std(a);
			b = b/std(b);
			x = xcorr(a,b); % positive peak means a lags b
			t = 1e-3*(1:length(x));
			t = t-mean(t);
			x = x/max(x);
			[~,loc] = max(x);
			peak_loc_xcorr_DA(j,i) = t(loc);

		end
	end
end


figure('outerposition',[0 0 600 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
clear l
l(1) = plot(mean_stim(:),peak_loc_xcorr(:)/1e-3,'k+');
l(2) = plot(mean_stim(:),peak_loc_xcorr_DA(:)/1e-3,'r+');

xlabel('Mean Stimulus (V)')
ylabel('Peak of xcorr (ms)')
legend(l,{'Data','DA Model'})

PrettyFig('plw=1.5;','lw=1.5;','fs=14;')

if being_published
	snapnow
	delete(gcf)
end

%   ########  ##     ## ##        ######  ########  ######             ##       ########     ###    
%   ##     ## ##     ## ##       ##    ## ##       ##    ##             ##      ##     ##   ## ##   
%   ##     ## ##     ## ##       ##       ##       ##                    ##     ##     ##  ##   ##  
%   ########  ##     ## ##        ######  ######    ######     #######    ##    ##     ## ##     ## 
%   ##        ##     ## ##             ## ##             ##              ##     ##     ## ######### 
%   ##        ##     ## ##       ##    ## ##       ##    ##             ##      ##     ## ##     ## 
%   ##         #######  ########  ######  ########  ######             ##       ########  ##     ## 

%%
% In this section we check if we can fit a DA model to the dose-response pulse experiment from Carlotta's paper. 

load('PulseData.mat')
d = struct;
for i = 1:length(PulseData)
	d(i).stimulus = mean2(PulseData(i).stim);
	d(i).response = mean2(PulseData(i).resp);
	d(i).response(3000:end) = [];
	d(i).stimulus(3000:end) = [];
	d(i).response(1:1000) = NaN;
	% d(i).response = [NaN(5,1); d(i).response];
	% d(i).stimulus = [zeros(5,1); d(i).stimulus];
end

% fit
p = FitModel2Data(@DAModelv2,d);


