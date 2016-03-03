% consolidateData.m
% consolidate data accepts a path, and merges all the .mat files there into one consolidate blob
% it makes many assumptions on how the data should be: it should be a Kontroller format
% and should have PID, voltage.
% this is not very general purpose, but hopefully will become more hardened to edge cases as time progresses
% 
% created by Srinivas Gorur-Shandilya at 11:09 , 06 July 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [PID, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes, sequence, all_spikes] = consolidateData(pathname,use_cache)

if ~nargin
	pathname = pwd;
end
if nargin < 2
	use_cache = false;
end

PID = [];
LFP = [];
fly = [];
spiketimes = [];
fA = [];
paradigm = [];
all_spikes = [];
orn = [];
AllControlParadigms = struct;
AllControlParadigms.Name = '';
AllControlParadigms.Outputs = [];
AllControlParadigms(1) = [];
paradigm_hashes = {};
sequence = []; % for each file. 

if ~strcmp(pathname(end),oss)
	pathname = [pathname oss];
end

allfiles = dir([pathname '*.mat']);
% remove the consolidated data from this
rm_this = [find(strcmp('cached_log.mat',{allfiles.name})) find(strcmp('consolidated_data.mat',{allfiles.name})) find(strcmp('cached.mat',{allfiles.name})) find(strcmp('template.mat',{allfiles.name}))];
if ~isempty(rm_this)
	allfiles(rm_this) = [];
end

% always be caching. 
hash = dataHash(allfiles);
if use_cache
	cached_data = [];
	try
		cached_data = load([pathname 'consolidated_data.mat'],'cached_data');
		cached_data = cached_data.cached_data;
	catch
	end
	if ~isempty(cached_data)
		PID = cached_data.PID;
		LFP = cached_data.LFP;
		fA = cached_data.fA;
		paradigm = cached_data.paradigm;
		AllControlParadigms = cached_data.AllControlParadigms;
		paradigm_hashes = cached_data.paradigm_hashes;
		orn = cached_data.orn;
		try
			all_spikes = cached_data.all_spikes;
		catch
		end
		try
			fly = cached_data.fly;
		catch
		end
		sequence = cached_data.sequence;
		return
	end
end

% find the longest trial in all the data. this will be the length of the combined data
disp('Determining longest data length...')
ll = 0;
for i = 1:length(allfiles)
	load(strcat(pathname,allfiles(i).name));
	for j = 1:length(data)
		if ~isempty(data(j).PID)
			ll = max([ll length(data(j).PID)]);
		end
	end
end
disp('Longest trial observed is:')
disp(ll)
ll = ceil(ll/10);
LFP = zeros(ll,0);
PID = zeros(ll,0);
fA = zeros(ll,0);
all_spikes = sparse(ll*10,0);

for i = 1:length(allfiles)
	clear spikes data metadata
	load(strcat(pathname,allfiles(i).name));
	disp(strcat(pathname,allfiles(i).name));
	for j = 1:length(data)
		clear this_PID this_LFP this_spikes this_hash this_paradigm 
		if ~isempty(data(j).PID)

			disp(['Paradigm : ' mat2str(j)]) 
			% figure out which control paradigm this is
			this_hash = dataHash(ControlParadigm(j));
			if isempty(find(strcmp(this_hash,paradigm_hashes)))
				AllControlParadigms(end+1) = orderfields(ControlParadigm(j));
				this_paradigm = length(AllControlParadigms);
				paradigm_hashes{end+1} = this_hash;
			else
				this_paradigm = find(strcmp(this_hash,paradigm_hashes));
			end

			disp('Loading PID and raw voltage...')
			this_PID = data(j).PID;
			this_LFP = data(j).voltage;
			this_fA = NaN*this_LFP(:,1:10:end)';

			disp(['There appear to be ' mat2str(width(this_PID)) ' trials for this paradigm.'])

			if exist('spikes','var')
				disp('spikes variable exists for this file.')
				if length(spikes) < j
					disp('No spike info for this paradigm.')
					this_spikes = sparse(size(this_LFP,1),size(this_LFP,2));
				else
					if length(spikes(j).A) > 10 && max(max(spikes(j).A)) > 0
						disp('Spikes exist for this paradigm, loading...')
						this_spikes = spikes(j).A;
						if width(this_spikes) < width(this_PID)
							disp('Spike width less than data width, padding...')
							this_spikes = sparse(size(this_LFP,1),size(this_LFP,2));
							for k = 1:width(this_spikes)
								try
									this_spikes(k,:) = spikes(j).A(k,:);
								catch
								end
							end
						end
					else
						disp('No spikes found!')
						% create dummy this_spikes
						this_spikes = sparse(size(this_LFP,1),size(this_LFP,2));
					end
				
					% censor trace use use_trace_fragment
					use_trace_fragment = [];
					try
						use_trace_fragment = spikes(j).use_trace_fragment;
					catch
					end

					if ~isempty(use_trace_fragment)
						disp('Censoring parts of trial based on annotation...')
						if width(use_trace_fragment) == width(this_PID)
							this_LFP(~logical(use_trace_fragment)) = NaN;
						else
							for k = 1:width(use_trace_fragment)
								this_LFP(k,~logical(use_trace_fragment(k,:))) = NaN;
							end
						end
					end
				end
			end

			% convert to firing rate
			this_fA = spiketimes2f(this_spikes,1e-4*(1:length(this_spikes)),1e-3,3e-2);

			rm_this = [];
			try
				rm_this = find(spikes(j).discard);
			catch
			end
			if ~isempty(rm_this)
				disp('Discarding some trials acc. to annotation:')
				this_PID(rm_this,:) = [];
				this_LFP(rm_this,:) = [];
				this_spikes(rm_this,:) = [];
				this_fA(:,rm_this) = [];

			else
				disp('No trials to discard')
			end

			disp('Downsampling data...')
			this_PID = this_PID(:,1:10:end)';
			this_LFP = this_LFP(:,1:10:end)';


			if length(this_LFP) < size(LFP,1) && length(this_LFP) > 1% this is to account for shorter paradigms, and pad them
				padding = NaN(size(LFP,1)-length(this_LFP),width(this_LFP));
				this_LFP = [this_LFP; padding];
				padding = NaN(size(PID,1)-length(this_PID),width(this_PID));
				this_PID = [this_PID; padding];
				padding = NaN(size(fA,1)-length(this_fA),width(this_fA));
				try
					this_fA = [this_fA; padding];
				catch
					this_fA = [this_fA(:); padding];
				end

				% pad spikes too
				padding = NaN(size(spikes,1)-length(this_spikes),width(this_spikes));
				this_spikes = [this_spikes; padding];
			end

			% check that PID and fA match
			if length(unique([width(this_PID) width(this_fA) width(this_LFP) width(this_spikes)])) > 1
				disp('Non matching PID and fA widths')
				keyboard
			end

			if length(unique([width(PID) width(fA) width(LFP) width(all_spikes)])) > 1
				disp('Non matching all_PID and fA widths')
				keyboard
			end

			% consolidate
			try
				LFP = [LFP this_LFP];
				PID = [PID this_PID];
				fA =  [fA  this_fA ];
				all_spikes = [all_spikes; this_spikes];
			catch er
				er
				keyboard
				
			end
			

			paradigm = [paradigm  this_paradigm*ones(1,width(this_PID))];
			orn = [orn  i*ones(1,width(this_PID))];

			% figure out the fly # from the file name
			try
				us = strfind(allfiles(i).name,'_'); % underscores
				this_fly = str2double(allfiles(i).name(strfind(allfiles(i).name,'_F')+2:us(find(us>strfind(allfiles(i).name,'_F')+2,1,'first'))-1));
				this_fly = 100*str2double(strrep(allfiles(i).name(1:us(3)-1),'_','')) + this_fly;
				fly = [fly; this_fly*ones(width(this_PID),1)];
			catch
				warning('Error in determining fly ID. Fly IDs may not match data.')
			end

			% also add the sequence 
			try
			this_sequence = find(timestamps(1,:)==j);
			if length(rm_this) 
				this_sequence(rm_this) = [];	
			end
			catch
				warning('Error in assembling sequence info.')
			end
			sequence = [sequence this_sequence];

		end
	end
end


% cache all of this
cached_data = struct;
cached_data.PID = PID;
cached_data.LFP = LFP;
cached_data.fA = fA;
cached_data.paradigm = paradigm;
cached_data.AllControlParadigms = AllControlParadigms;
cached_data.paradigm_hashes = paradigm_hashes;
cached_data.orn = orn;
cached_data.fly = fly;
cached_data.all_spikes  = all_spikes;
cached_data.sequence = sequence;
save([pathname 'consolidated_data.mat'],'cached_data');

