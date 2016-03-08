% raw2ORNData.m
% converts a folder with kontroller-generated files into ORNData class objects, one for each file
% uses consolidateData for some intermediate processing 
% 
% created by Srinivas Gorur-Shandilya at 3:20 , 02 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [od] = raw2ORNData(pathname,varargin)


% options and defaults
options.use_led = false;
options.filter_LFP = false; 
options.led_power_func = [];

assert(ischar(pathname),'First argument should be a string specifying the path to the raw data files')
assert(exist(pathname,'dir')==7,'Path not found')
assert(exist('consolidateData','file')==2,'consolidateData.m not found');

% validate and accept options
if iseven(length(varargin))
	for ii = 1:2:length(varargin)-1
	temp = varargin{ii};
    if ischar(temp)
    	if ~any(find(strcmp(temp,fieldnames(options))))
    		disp(['Unknown option: ' temp])
    		disp('The allowed options are:')
    		disp(fieldnames(options))
    		error('UNKNOWN OPTION')
    	else
    		options = setfield(options,temp,varargin{ii+1});
    	end
    end
end
else
	error('Inputs need to be name value pairs')
end

[stimulus, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes, sequence, all_spikes] = consolidateData(pathname,true);

% use the LED if necessary
if options.use_led
	stimulus = NaN*fA;
	for i = 1:length(paradigm)
		temp = AllControlParadigms(paradigm(i)).Name;
		stimulus(:,i) = options.led_power_func(AllControlParadigms(paradigm(i)).Outputs(1,1:10:end));
	end
end

% get a list of all the file names
if ~strcmp(pathname(end),oss)
	pathname = [pathname oss];
end

allfiles = dir([pathname '*.mat']);
% remove the consolidated data from this
rm_this = [find(strcmp('cached_log.mat',{allfiles.name})) find(strcmp('consolidated_data.mat',{allfiles.name})) find(strcmp('cached.mat',{allfiles.name})) find(strcmp('template.mat',{allfiles.name}))];
if ~isempty(rm_this)
	allfiles(rm_this) = [];
end

% remove some junk
rm_this = max(fA) == 0;
fA(:,rm_this) = [];
stimulus(:,rm_this) = [];
LFP(:,rm_this) = [];
orn(rm_this) = [];
paradigm(rm_this) = [];

% do we have to filter the LFP?
if options.filter_LFP
	filtered_LFP = LFP;
	for i = 1:width(LFP)
		filtered_LFP(:,i) = filtered_LFP(:,i) - fastFiltFilt(ones(1e4,1),1e4,filtered_LFP(:,i));
		filtered_LFP(:,i) = filtered_LFP(:,i)*10; % to get the units right, now in mV
	end
else
	filtered_LFP = LFP*10;
end

% organise by paradigm and by orn
for i = max(orn):-1:1
	for j = max(paradigm):-1:1
		od(i,j) = ORNData;
		od(i,j).original_name = allfiles(i).name;
		od(i,j).firing_rate = fA(:,orn==i & paradigm == j);
		od(i,j).stimulus = stimulus(:,orn==i & paradigm == j);
		od(i,j).LFP = filtered_LFP(:,orn==i & paradigm == j);
	end
end
