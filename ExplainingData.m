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


% Fit DA Model to Mean Shifted Gaussians 


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
end


