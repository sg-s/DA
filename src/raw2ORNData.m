% raw2ORNData.m
% converts a folder with kontroller-generated files into ORNData class objects, one for each file
% uses consolidateData for some intermediate processing 
% 
% created by Srinivas Gorur-Shandilya at 3:20 , 02 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [od] = raw2ORNData(pathname,light_power_fit)

if ~nargin
	pathname = pwd;
end
if nargin < 2
	use_led = false;
else
	use_led = true;
end

assert(exist('consolidateData','file')==2,'consolidateData.m not found');

[stimulus, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes, sequence, all_spikes] = consolidateData(pathname,true);

% use the LED if necessary
if use_led
	stimulus = NaN*fA;
	for i = 1:length(paradigm)
		temp = AllControlParadigms(paradigm(i)).Name;
		stimulus(:,i) = light_power_fit(AllControlParadigms(paradigm(i)).Outputs(1,1:10:end));
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

for i = 1:max(orn)
	od(i) = ORNData;
	od(i).original_name = allfiles(i).name;
	od(i).paradigm = paradigm(orn==i);
	od(i).firing_rate = fA(:,orn==i);
	od(i).stimulus = stimulus(:,orn==i);
	od(i).LFP = LFP(:,orn==i);
end
