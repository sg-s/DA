% find fudge factor

fudge_factor_range = [0.5:0.1:1];

% first, use manual control to make sure you have depleted the headspace
% nicely. 
alldata = struct;

for i = 1:length(fudge_factor_range)
    fudge_factor = fudge_factor_range(i);
    
    % make control paradigm with this
    ControlParadigm = makeVarianceSwitching(fudge_factor);
    
    % now run the experiment
    alldata(i).data = kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',[2 2],'w',10000);
    
end