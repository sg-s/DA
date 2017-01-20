% get all files

disp('Getting full file list...')
if ~exist('allfiles','var')
	allfiles = getAllFiles(pwd);
end
disp('Done!')

% for each file, convert it to TDMS and save it in the kontroller format
for i = 1:length(allfiles)
	[~,~,ext] = fileparts(allfiles{i});
	if strcmp(ext,'.mat') & ~any(strfind(allfiles{i},'/.'))
		clear v

			v = whos('-file',allfiles{i});
		
			if any(strcmp({v.name},'ConvertedData'))
				disp('Converting:')
				disp(allfiles{i})
				clear data ConvertedData ControlParadigm data metadata OutputChannelNames SamplingRate timestamps 
				load(allfiles{i})

				metadata.info = 'converted from Carlotta data';
				timestamps = [];
				OutputChannelNames = {'dummy'};

				% figure out the sampling rate
				if length(ConvertedData.Data.Root.Property) == 1
					if strcmp(ConvertedData.Data.Root.Property.Name,'sampling rate')
						SamplingRate = ConvertedData.Data.Root.Property.Value;
					else 
						error('unknown property, cannot determine sampling rate') 	
					end
				else
					error('More than one property')
				end

				for j = 1:length(ConvertedData.Data.MeasuredData)
					% figure out which trial, channel it is
					t = ConvertedData.Data.MeasuredData(j).Name;
					if any(strfind(t,'trial')) & any(strfind(t,'/'))
						first_space = (strfind(t,' ')); first_space = first_space(1);
						slash_loc = strfind(t,'/');
						trial_no = str2double(t(first_space:slash_loc-1));
						this_channel = strtrim(t(slash_loc+1:end));
						if strcmp(this_channel,'response')
							this_channel = 'voltage';
						end
						data(trial_no).(this_channel) = ConvertedData.Data.MeasuredData(j).Data';
						ControlParadigm(trial_no).Outputs = 0*ConvertedData.Data.MeasuredData(j).Data';
						ControlParadigm(trial_no).Name = 'dummy';
					end
				end



				% save
				new_name = [allfiles{i}(1:end-3) 'kontroller'];
				save(new_name,'data','ControlParadigm','timestamps','SamplingRate','OutputChannelNames','metadata');

				% delete the old one
				delete(allfiles{i})
			end


	end
end