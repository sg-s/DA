% compressMMimages.m
% batch converts an image sequence into a mat file
% created by Srinivas Gorur-Shandilya at 9:21 , 03 December 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function compressMMimages(source_path,destination_path)

if ~nargin
	source_path = 'Z:\srinivas_gs\calcium-images\raw-images\';
	destination_path = 'Z:\srinivas_gs\calcium-images\compressed-images\';
end

if ~strcmp(destination_path(end),oss)
	destination_path = [destination_path oss];
end

all_sub_folders = getAllSubFolders(source_path);



for i = 1:length(all_sub_folders)
	allfiles = dir([all_sub_folders{i} '*.tif']);
	if length(allfiles)
		% read the metadata
		txt = fileread([all_sub_folders{i} 'metadata.txt']);
		hash = dataHash(txt);

		% look for a filename with this hash in the destination_path
		compressed_files = dir(destination_path);
		if any(find(strcmp(['video_' hash '.mat'],{compressed_files.name})))
			disp('Already compressed. Skipping...')
		else
			disp(hash)
			disp('New data. Compressing...')

			imgSeq2mat(all_sub_folders{i});

			% now move the compressed mat file to the destination
			movefile([all_sub_folders{i} 'video_' hash '.mat'],[destination_path 'video_' hash '.mat'])

		end

		
	end
end