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

function [PID, LFP, fA, paradigm, orn, AllControlParadigms, paradigm_hashes] = consolidateData(pathname,use_cache)


PID = [];
LFP = [];
fA = [];
paradigm = [];
orn = [];
AllControlParadigms = struct;
AllControlParadigms.Name = '';
AllControlParadigms.Outputs = [];
AllControlParadigms(1) = [];
paradigm_hashes = {};

allfiles = dir([pathname '*.mat']);

% always be caching. 
hash = DataHash(allfiles);
if use_cache
	cached_data = cache(hash);
	if ~isempty(cached_data)
		PID = cached_data.PID;
		LFP = cached_data.LFP;
		fA = cached_data.fA;
		paradigm = cached_data.paradigm;
		AllControlParadigms = cached_data.AllControlParadigms;
		paradigm_hashes = cached_data.paradigm_hashes;
		orn = cached_data.orn;
		return
	end
end

for i = 1:length(allfiles)
	load(strcat(pathname,allfiles(i).name));
	for j = 1:length(data)
		if ~isempty(data(j).PID)
			clear this_paradigm 
			% figure out which control paradigm this is
			this_hash = DataHash(ControlParadigm(j));
			if isempty(find(strcmp(this_hash,paradigm_hashes)))
				AllControlParadigms(end+1) = ControlParadigm(j);
				this_paradigm = length(AllControlParadigms);
				paradigm_hashes{end+1} = this_hash;
			else
				this_paradigm = find(strcmp(this_hash,paradigm_hashes));
			end


			this_PID = data(j).PID;
			this_LFP = data(j).voltage;
			this_fA = NaN*this_LFP(:,1:10:end)';

			if length(spikes) < j
			else
				if ~isempty(spikes(j).A)
					this_fA = spiketimes2f(spikes(j).A,1e-4*(1:length(spikes(j).A)),1e-3,3e-2);
				end
				

				% censor trace use use_trace_fragment
				use_trace_fragment = [];
				try
					use_trace_fragment = spikes(j).use_trace_fragment;
				catch
				end

				if ~isempty(use_trace_fragment)
					if width(use_trace_fragment) == width(this_PID)
						this_LFP(~logical(use_trace_fragment)) = NaN;
					else
						for k = 1:width(use_trace_fragment)
							this_LFP(k,~logical(use_trace_fragment(k,:))) = NaN;
						end
					end
				end


				rm_this = [];
				try
					rm_this = find(spikes(j).discard);
				end
				if ~isempty(rm_this)
					this_PID(rm_this,:) = [];
					this_LFP(rm_this,:) = [];
					try
						this_fA(:,rm_this) = [];
					catch
					end
				else

				end

			
			end

			this_PID = this_PID(:,1:10:end)';
			this_LFP = this_LFP(:,1:10:end)';


			if width(this_fA) ~= width(this_LFP)
				this_fA = [this_fA NaN(length(this_fA),width(this_LFP) - width(this_fA))];
			end

			% consolidate
			LFP = [LFP this_LFP];
			PID = [PID this_PID];
			fA =  [fA  this_fA ];
			

			paradigm = [paradigm  this_paradigm*ones(1,width(this_PID))];
			orn = [orn  i*ones(1,width(this_PID))];

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
cache(hash,[]);
cache(hash,cached_data);


