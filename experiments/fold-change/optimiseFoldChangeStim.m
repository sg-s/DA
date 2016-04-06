% this script tunes the MFCs to make fold change stimuli
% target stimuli:
% 0.1 -> 0.2 V
% 0.2 -> 0.4 V
% 0.4 -> 0.8 V
% 0.8 -> 1.6 V

MFC200_setpoints = [1 1 1 1];
MFC500_setpoints = [1 1 1 1];

MFC200_max = 5;
MFC500_max = 5;

MFC200_min = .1;
MFC500_min = .1;

nsteps = 50;

% make the figure
f = figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
for i = 6:-1:1
	ax(i) = subplot(2,3,i); hold on
end

for i = 1:nsteps
	suptitle(['Iteration #' oval(i)])

	% make the stimuli
	ControlParadigm = makeFoldChangeStim(MFC200_setpoints,MFC500_setpoints);

	% run the trial
	data = kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',[1 2 3 4],'w',1e4);

	% for each paradigm, get a score of how bad you are doing and adjust your setpoints
	for j = 1:length(data)
		
	end


end
