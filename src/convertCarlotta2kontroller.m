% get all files

allfiles = getAllFiles(pwd);

% for each file, convert it to TDMS and save it in the kontroller format
for i = 1:length(allfiles)
	[~,~,ext] = fileparts(allfiles{i});
	if strcmp(ext,'.tdms')
		disp(allfiles{i})
		try
			convertTDMS(1,allfiles{i});
			delete(allfiles{i})
		catch
		end
	end
end