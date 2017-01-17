% this script searches all of carlotta's raw data for files where the PID and the firing rate are co-recorded 

disp('Getting list of all files to scan...')
allfiles = getAllFiles;

haz_stimulus = false(length(allfiles),1);
haz_response = false(length(allfiles),1);
PID_snippets = NaN(2e3,length(allfiles));

figure, hold on

for i = 1:length(allfiles)
	%textbar(i,length(allfiles))

	[~,~,ext] = fileparts(allfiles{i});
	if strcmp(ext,'.tdms')
		[ConvertedData,~,ChanNames]=convertTDMS(0,allfiles{i});
		ChanNames = ChanNames{1};
		if any(cell2mat((cellfun(@(x) strfind(x,'PID'),ChanNames,'UniformOutput',false))))
			haz_stimulus(i) = true;
		end

		if any(cell2mat((cellfun(@(x) strfind(x,'response'),ChanNames,'UniformOutput',false))))
			haz_response(i) = true;
		end

		% grab all the PID
		this_PID = NaN(2e3,length(ConvertedData.Data.MeasuredData));
		for j = 1:length(ConvertedData.Data.MeasuredData)
			if any(strfind(ConvertedData.Data.MeasuredData(j).Name,'PID'))
				[~,loc] = max(ConvertedData.Data.MeasuredData(j).Data);
				a = loc - 1e4;
				z = loc + 1e4 - 1;
				try
					temp = ConvertedData.Data.MeasuredData(j).Data(a:z);
					this_PID(:,j) = temp(1:10:end);
					
					idx = find(max(this_PID) == max(max(this_PID)),1,'first');
					PID_snippets(:,j) = this_PID(:,idx);
					plot(PID_snippets(:,j))
					drawnow
				catch
				end
			end
		end


	end

end

