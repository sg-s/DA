% this is carlotta's function
function [timesr PIDs sr msr sdr] = bin_traces(PID, stA, deltat, win, sliding)

% bins PID single traces
[num_dil, trials, L] = size(PID);
PIDs =  [];
for dil = 1:num_dil
    for t = 1:trials
        [timesr PIDs(dil,t,:)] = sliding_average(squeeze(PID(dil, t,:)), deltat, win, sliding);
        PIDs(dil,t,:) = PIDs(dil,t,:)-mean(PIDs(dil,t,find(timesr>1,1):find(timesr>4,1)));
        PIDs(dil,t,:) = PIDs(dil,t,:)/max(PIDs(dil,t,:));
    end
end

if sum(any(stA))>0
    % bins spikes in single traces
    [num_dil, trials, maxtime] = size(stA);
    sr = [];
    maxtime = L;
    for dil = 1:num_dil
        for t = 1:trials
            times = squeeze(stA(dil,t,:));
            if any(times)
                times(times==0)=[];
                spk = spiketime2spk(times',maxtime);
                [timesr sr(dil,t,:)] = spike_rate(spk', deltat, win, sliding);
            end
        end
    end
    
    % calculate mean firing rate
    msr = zeros(num_dil, length(sr));
    sdr = zeros(num_dil, length(sr));
    for dil=1:num_dil
        SR = squeeze(sr(dil,:,:));
        jj=0;
        good = [];
        for t=1:size(sr,2)
            if any(SR(t,:))
                jj=jj+1;
                good(jj) = t;
            else
                display('found empty trace')
            end
        end
        msr(dil, :) = mean(SR(good,:));
        sdr(dil, :) = std(SR(good,:));
        
    end
else
   msr = 0;
   sdr = 0;
   sr = 0;
end



end

