% consolidateData2.m
% consolidate data accepts a path, and merges all the .mat files there into one consolidate blob
% it makes many assumptions on how the data should be: it should be a Kontroller format
% and should have PID, voltage.
% this is not very general purpose, but hopefully will become more hardened to edge cases as time progresses
% 
% created by Srinivas Gorur-Shandilya at 11:09 , 06 July 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function consolidated_data = consolidateDataPairedPulses(pathname,varargin)

if ~nargin
	pathname = pwd;
end

options.use_cache = true;


if nargout && ~nargin 
	varargout{1} = options;
    return
end

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
    		options.(temp) = varargin{ii+1};
    	end
    end
end
elseif isstruct(varargin{1})
	% should be OK...
	options = varargin{1};
else
	error('Inputs need to be name value pairs')
end

PID = [];
fly = [];
paradigm = [];
orn = [];
fA = [];


% paired pulse specific
pulse_seperation = [];
background_id = [];

if ~strcmp(pathname(end),oss)
	pathname = [pathname oss];
end

allfiles = dir([pathname '*.mat']);
% remove the consolidated data from this
rm_this = [find(strcmp('cached_log.mat',{allfiles.name})) find(strcmp('consolidated_data.mat',{allfiles.name})) find(strcmp('cached.mat',{allfiles.name})) find(strcmp('template.mat',{allfiles.name}))];
if ~isempty(rm_this)
	allfiles(rm_this) = [];
end

% remove stupid ass shit created by fucking Dropbox
rm_this = cell2mat(cellfun(@(x) any(strfind(x,'._')), {allfiles.name},'UniformOutput',false));
allfiles(rm_this) = [];

% always be caching. 
if options.use_cache
	consolidated_data = [];
	try
		consolidated_data = load([pathname 'consolidated_data.mat'],'consolidated_data');
		consolidated_data = consolidated_data.consolidated_data;
	catch
	end
	if ~isempty(consolidated_data)
		return
	end
end

% find the longest trial in all the data. this will be the length of the combined data
disp('Determining longest data length...')
ll = 0;
for i = 1:length(allfiles)
	load(strcat(pathname,allfiles(i).name));
	for j = 1:length(data)
		if ~isempty(data(j).voltage)
			ll = max([ll length(data(j).voltage)]);
		end
	end
end

disp('Longest trial observed is:')
disp(ll)
ll = ceil(ll/10); % because we want everything at a 1ms timestep

% make placeholders
LFP = zeros(ll,0);
PID = zeros(ll,0);
fA = zeros(ll,0);


for i = 1:length(allfiles)
	clear spikes data metadata
	load(strcat(pathname,allfiles(i).name));
	disp(strcat(pathname,allfiles(i).name));

	assert(length(spikes) == length(data),'Data and spikes are of different sizes')

	% concat everythign
	this_PID = vertcat(data.PID);
	this_LFP = vertcat(data.voltage);

	% rotate, and subsample
	this_LFP = this_LFP(:,1:10:end)';
	this_PID = this_PID(:,1:10:end)';


	% calculate the A firing rates
	this_fA = zeros(ll,0);
	for j = 1:length(spikes)
		temp = spiketimes2f(spikes(j).A,1e-4*(1:length(spikes(j).A)),1e-3,3e-2);
		this_fA = [this_fA temp];
	end

	assert(size(this_fA,2) == size(this_PID,2),'unequal stimulus and response sizes!')	

	% get fly ID
	us = strfind(allfiles(i).name,'_'); % underscores
	this_fly = str2double(allfiles(i).name(strfind(allfiles(i).name,'_F')+2:us(find(us>strfind(allfiles(i).name,'_F')+2,1,'first'))-1));
	this_fly = 100*str2double(strrep(allfiles(i).name(1:us(3)-1),'_','')) + this_fly;
	this_fly = this_fly*ones(size(this_LFP,2),1);

	% get sensillum ID
	this_orn = i*ones(size(this_LFP,2),1);

	% measure inter-pulse distance for each row here 
	these_paradigms = NaN*ones(size(this_LFP,2),1);
	a = 1;
	for j = 1:length(data)
		s = size(data(j).voltage,1);
		z = a + s - 1;
		these_paradigms(a:z) = j;
		a = z + 1;
	end

	this_pulse_sep = NaN*these_paradigms;
	for j = 1:length(these_paradigms)
		this_paradigm = these_paradigms(j);
		[ons,offs]=computeOnsOffs(ControlParadigm(this_paradigm).Outputs(4,:));
		this_pulse_sep(j) = 1e-4*(ons(2) - offs(1));
	end


	% merge
	LFP = [LFP this_LFP];
	PID = [PID this_PID];
	fA = [fA  this_fA];
	pulse_seperation = [pulse_seperation; this_pulse_sep(:)];
	orn = [orn; this_orn(:)];
	fly = [fly; this_fly(:)];


end


% cache all of this
consolidated_data = struct;
consolidated_data.PID = PID;
consolidated_data.fA = fA;
consolidated_data.LFP = LFP;
consolidated_data.orn = orn;
consolidated_data.fly = fly;

consolidated_data.pulse_seperation = pulse_seperation;
consolidated_data.background_id = background_id;
save([pathname 'consolidated_data.mat'],'consolidated_data');

