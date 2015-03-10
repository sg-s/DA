% ConvertCarlottaFig3
% converts carlotta's raw data in Fig 3 (a,b,c) of her paper
% to a Kontroller-friendly form
% removes weird data, cleans up various things
% 
% created by Srinivas Gorur-Shandilya at 3:58 , 17 February 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

filepath = ('/local-data/DA-paper/carlotta/fig3/data_2011_10_04_2ac_dr_N12.mat');

o.Input = 'File';
hash = DataHash(filepath,o);
if ~strcmp('f550ece607be1583c622aad49e93d060',hash)
	error('Hash failure with data. Are you sure you''re analysing the right file?')
end

clear data ControlParadigm
load(filepath)

for i = 1:length(dilutions)
	data(i).PID = squeeze(PID(i,:,:));
	data(i).Voltage = squeeze(ORN(i,:,:));
	ControlParadigm(i).Name = oval(log10(dilutions(i)));
	ControlParadigm(i).Outputs = zeros(1,length(data(i).PID));
end


SamplingRate = 1e4;
OutputChannelNames = {};
save('/local-data/DA-paper/carlotta/fig3/abc.mat','ControlParadigm','data','SamplingRate','OutputChannelNames')

