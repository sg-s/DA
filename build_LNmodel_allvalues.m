% this does NOT build the LN model--it makes the prediction based on an already calcualted filter
% all the useful info is in model_values, which I explicitly pass back
function [fSTA, fLN, model_values] = build_LNmodel_allvalues(time2, stim, f2, fsd2, trials, memory, deltat, fitfun, figon, a, D)

L = length(stim);

s = zeros(L-memory, memory+1);
for i=1:L-memory
    s(i,:) = stim(memory+i:-1:i);
end
fSTA = s*a;
    
if fitfun==1 % fit hill function
    fLN = D(1)*(1./(1+(D(2)./(fSTA-D(4))).^D(3)));
    actualfp = fLN; % we will return this with model values...
    fLN(fSTA<D(4))=0; % need this to prevent complex numbers coming out from the hill function
else % fit line/rectifier
    fLN = D(1)*fSTA + D(2);
    zeroresp = min(f2);
    % Here Carlotta is throing away all predicted firing rates lower than the minimum firing rate observered in the data. cheating...
    actualfp = fLN; % we will return this with model values...
    fLN = fLN.*(fLN>zeroresp) + zeroresp*(fLN<zeroresp);
end
%         Pr = mean(abs(f2-fLN).^2);
%         Pn = mean(fsd2.^2);
%         NR = sqrt(Pr/Pn);

if figon==1
    figure;
    
    subplot(2,2,1)
    plot((1:memory+1)*deltat, a, 'r', 'linewidth', 2)
    
    
    subplot(2,2,2); hold on
%    plot(fSTA, f2, '.', 'color', [.85 .85 .85])
    %             plot(bfSTA(1:nbins)+mfSTA, N, 'k', 'linewidth', 2)
    if fitfun==1
        ffitted = @(x) D(1)*(1/(1+(D(2)/(x-D(4)))^D(3)));
    else 
        ffitted = @(x) D(1)*x + D(2); 
    end
    

    fplot(ffitted, [min(fSTA) max(fSTA)], 'r');
    
    subplot(2,2,3:4); hold on
    errorbar_cazzo(time2, f2, fsd2/sqrt(trials))

    plot(time2, fLN, 'r')
end

model_values.filtertime = (1:memory+1)*deltat;
model_values.K = a; % filter
if fitfun==1
        ffitted = @(x) D(1)*(1/(1+(D(2)/(x-D(4)))^D(3)));
    else 
        ffitted = @(x) D(1)*x + D(2); 
    end
model_values.NLfunction = ffitted;
model_values.D = D;
model_values.NLfunctionRange = [min(fSTA) max(fSTA)];
model_values.time = time2;
model_values.f = f2;
model_values.fs = fsd2/sqrt(trials);
model_values.fp = fLN;
model_values.actualfp = actualfp; % actual prediction, including firing rates < 0
    
    
end

