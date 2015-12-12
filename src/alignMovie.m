function [output, xShifts, yShifts, scoreChosen, scoreBestNoLoss, positionBestNoLoss] = AlignMovie(input, target)

inputSingle = single(input);
clear input;
%inputSingle = inputSingle;
[sizeX, sizeY, numTimepoints] = size(inputSingle);

if nargin < 2
    subArray = inputSingle(:,:,10:50);
    target = mean(subArray,3);
end
threshold = prctile(reshape(target,1,[]),65);
targetBinary = single(target > threshold);
targetGPU = (targetBinary);

% gaussianFilterZ = gausswin(3,2);
% gaussianFilterZ = gaussianFilterZ/(sum(gaussianFilterZ));
% gaussianFilterXY = gausswin(10,2);
% gaussianFilterXY = gaussianFilterXY/(sum(gaussianFilterXY));
% gpuInputConv = convn(gpuInput,(permute(gaussianFilterZ,[3 2 1])),'same');
% gpuInputConv = convn(gpuInputConv,(gaussianFilterXY),'same');
% gpuInputConv = convn(gpuInputConv,(gaussianFilterXY'),'same');
chunkSize = 500;
numChunks = ceil(numTimepoints/chunkSize);
gaussFilterSigma = [2.5,2.5,2.5];
gaussFilterSize = [11,11,11];
zOverlap = ceil((gaussFilterSize(3)-1)/2);

crossCorrs(2*sizeX-1,2*sizeY-1,numTimepoints) = single(0);
fprintf('%4d%% done', 0)
for chunk = 1:numChunks
    %chunkSize is the number of frames that will be correctly calculated in
    %each iteration. validChunkStart and End are inclusive, MATLAB indexing
    %style.
    validChunkStart = (chunk-1)*chunkSize+1;
    validChunkEnd = chunk*chunkSize;
    if validChunkEnd > numTimepoints
        validChunkEnd = numTimepoints;
    end
    
    %In order to calculate these values, we need to grab frames zOverlap
    %before and after the valid chunk
    inputChunkStart = validChunkStart - zOverlap;
    if inputChunkStart < 1
        inputChunkStart = 1;
    end
    inputChunkEnd   = validChunkEnd   + zOverlap;
    if inputChunkEnd > numTimepoints
        inputChunkEnd = numTimepoints;
    end
    
    %When outputting, we need to get rid of these excess overlap values
    numCalculatedTimepoints = (1+inputChunkEnd-inputChunkStart);
    if chunk == 1
        outputChunkStart = 1;
    else
        outputChunkStart = 1 + zOverlap;
    end
    if chunk == numChunks
        outputChunkEnd = numCalculatedTimepoints;
    else
        outputChunkEnd = numCalculatedTimepoints - zOverlap;
    end
    

    gpuInput = (inputSingle(:,:,inputChunkStart:inputChunkEnd));
    gpuInputConv = imgaussfilt3(gpuInput,gaussFilterSigma,'FilterSize',gaussFilterSize);
    clear gpuInput;

    blurred = gather(gpuInputConv);
    clear gpuInputConv;
    binary(sizeX,sizeY,numCalculatedTimepoints) = single(0);
    for t = 1:numCalculatedTimepoints
        threshold = prctile(reshape(blurred(:,:,t),1,[]),65);
        binary(:,:,t) = single(blurred(:,:,t) > threshold);
    end
    binaryGPU = (binary);
    clear binary;

    %Allocate the correct amount of space on GPU by grabbing from CPU
    crossCorrsGPU = (crossCorrs(:,:,inputChunkStart:inputChunkEnd));


    for t = 1:numCalculatedTimepoints
        crossCorrsGPU(:,:,t) = normxcorr2(binaryGPU(:,:,t),targetGPU);
        if (~mod(t,100))
            chunksDone = round(100*(t+(chunk-1)*chunkSize)/(numTimepoints));
            fprintf('\b\b\b\b\b\b\b\b\b\b%4d%% done', chunksDone)
%             disp(t+(chunk-1)*chunkSize)
        end
    end
    clear binaryGPU;

    crossCorrs(:,:,validChunkStart:validChunkEnd) = gather(crossCorrsGPU(:,:,outputChunkStart:outputChunkEnd));
    clear crossCorrsGPU;
    
end
fprintf('\n');

[Y, X] = meshgrid(1:(2*sizeY-1),1:(2*sizeX-1));
centeredY = Y-sizeY;
centeredX = X-sizeX;
lengthScale = sqrt((.5*sizeX)^2 + (.5*sizeY)^2);
distanceFromCenter = sqrt((centeredX).^2 + (centeredY).^2)/lengthScale;
centerDistancePenalty = (-nthroot(2*(distanceFromCenter-0.5),3)+1)/2;% -log(log(-distanceFromCenter+exp(1))).^2+1;

scoreChosen(numTimepoints, 1) = 0;
scoreBestNoLoss(numTimepoints, 1) = 0;
positionBestNoLoss(numTimepoints, 2) = 0;
xShifts(numTimepoints, 1) = 0;
yShifts(numTimepoints, 1) = 0;

fprintf('Finding optimal alignment\n');
fprintf('%4d%% done', 0)
for t = 1:numTimepoints
    if (t > 1)
        previousXYPos = xyPos(:,t-1);
    else
        previousXYPos = [sizeX sizeY];
    end
    distanceFromPrevious = sqrt((X-previousXYPos(1)).^2 + (Y-previousXYPos(2)).^2)/lengthScale;
    previousDistancePenalty = distanceFromPrevious<lengthScale;%(-nthroot(2*(distanceFromPrevious-0.5),3)+1)/2;% -log(log(-distanceFromPrevious+exp(1))).^2+1;

    [scoreBestNoLoss(t), indBest] = max(reshape(crossCorrs(:, :, t), 1, []));
    [positionBestNoLoss(t, 1) positionBestNoLoss(t, 2)] = ind2sub(size(crossCorrs(:, :, t)), indBest);
    positionBestNoLoss(t, :) = positionBestNoLoss(t, :) - [sizeX sizeY];
    
    score = crossCorrs(:,:,t) .* centerDistancePenalty .* previousDistancePenalty;
%     score = crossCorrs(:,:,t);
    maxScore = max(max(score));
    scoreChosen(t) = maxScore;
    [xPos, yPos] = find(score == maxScore,1);
    xyPos(:,t) = [xPos; yPos];
    
    xyShifts = xyPos(:,t) - [sizeX; sizeY];
    xShifts(t) = xyShifts(1);
    yShifts(t) = xyShifts(2);
    if (~mod(t,100))
            chunksDone = round(100*t/(numTimepoints));
            fprintf('\b\b\b\b\b\b\b\b\b\b%4d%% done', chunksDone)
%             disp(t+(chunk-1)*chunkSize)
    end
end
fprintf('\n');

clear crossCorrs;

output(sizeX,sizeY,numTimepoints) = 0;
for t = 1:numTimepoints
    output(:,:,t) = ShiftMatrix(inputSingle(:,:,t),1,[yShifts(t),xShifts(t)],0);
end