function result = runningVote(filename, withFilter )

warning('off', 'Images:initSize:adjustingMag');
% clear variables;
% close all;

correctEncoding = decodeFilename(filename);
% fprintf('Original: %s\n', sprintf('%d', correctEncoding ) );

%set parameters for acceleration
codeLength = 24;
numSamplesPerBar = 500;
runningVoteThreshold = 64;
confidenceThreshold = 0.8;
searchStencilWidth = 181; 
directionStencilWidth = 21;
nHood = strel('diamond',3);
vectorsAlignedThresh = 0.35;

stencil = zeros(searchStencilWidth,searchStencilWidth);
stencil( (searchStencilWidth+1)/2, (searchStencilWidth+1)/2 ) = 1;
stencil = bwdist(stencil);

smallStencil = zeros(directionStencilWidth,directionStencilWidth);
smallStencil( (directionStencilWidth+1)/2, (directionStencilWidth+1)/2) = 1;
smallStencil = bwdist(smallStencil);

%% labelPixels
[labeled, nred, nblue] = labelRBPix( filename, withFilter );
tmpLabeled = labeled;
tmpLabeled( labeled == 0 ) = -(nred+nblue);
% figure; imshow(tmpLabeled, []);

%% Create Graph of Neighbors and Normal-Directions
numlabels = nred + nblue;
floodGraph = cell(numlabels, 5);
for label = 1:numlabels
    % Graph of Neighbors
    currLabelPix = labeled == label;
    neighbors = imdilate(currLabelPix, nHood) & ~currLabelPix;
    
    neighborLabels = unique( labeled(neighbors) );
    neighborLabels( neighborLabels == 0 ) = []; % remove boundary
    floodGraph{label, 1} = neighborLabels;
    
    % Directions between neighbors
    labeledNeighbors = maskImage( labeled, neighbors );    
    directionsFromLabelToNeighbor = zeros( length(neighborLabels), 2 );
    for n = 1:length(neighborLabels)
        neighbor = neighborLabels(n);
        
        [rows, cols] = find( labeledNeighbors == neighbor );
        numSamples = min( numSamplesPerBar, length(rows));
        index = randsample(1:length(rows), numSamples);
        rows = rows(index,:);
        cols = cols(index,:);

        dirs = zeros(numSamples, 2);
        validSamples = 0;
        for s = 1:numSamples
            center = [rows(s) cols(s)];
            pixelFound = findWithStencil( labeled, center, label, smallStencil);
            if pixelFound(1) ~= -1
               validSamples = validSamples+1;
               dir = center - pixelFound;
               dirs(validSamples,:) = dir / norm(dir);
            end
        end
        dirs = dirs(1:validSamples,:);
        dirFromLabelToNeighbor = mean(dirs,1);
        directionsFromLabelToNeighbor(n,:) = dirFromLabelToNeighbor;
    end
    floodGraph{label, 3} = directionsFromLabelToNeighbor;
end

%% Compute pairwise ratios via local points
startLabels = [];
distPairMap = containers.Map('KeyType','char','ValueType','double');
for label = 1:numlabels
    neighborLabels = floodGraph{label,1};
    
    for n = 1:length(neighborLabels)
        neighbor = neighborLabels(n);
        keyPair = sprintf('%d_%d', sort([label, neighbor]));
        if ~isKey(distPairMap, keyPair)

            % ignore ratio if either region is on border
            neighborNeighborLabels = floodGraph{neighbor,1};
            if length( neighborLabels ) < 2 || length( neighborNeighborLabels ) < 2
                distPairMap(keyPair) = NaN;
                continue;
            end            
            
            % get shared points between region-pair
            currLabelPix = labeled == label;
            neighbors = imdilate(currLabelPix, nHood) & ~currLabelPix;
            labelBorder = neighbors & bwperim( labeled == neighbor );
            currLabelPix = labeled == neighbor;
            neighborBorder = imdilate(currLabelPix, nHood) & ~currLabelPix & bwperim( labeled == label );
            sharedBorder = labelBorder | neighborBorder;

            [rows, cols] = find( sharedBorder == 1 );
            numSamples = min( numSamplesPerBar, length(rows));
            index = randsample(1:length(rows), numSamples);
            rows = rows(index,:);
            cols = cols(index,:);
            localRatioPairs = zeros(numSamples, 2);

            % compute lengths at points for Label to other neighbors
            neighborsBorderMasked = maskImage( labeled, neighbors );
            neighborsBorderMasked( neighborsBorderMasked == neighbor ) = 0;
            dirFromLabelToNeighbor = floodGraph{label,3}(n,:);
            
            % ignore neighbors if would cause loop            
            for v = 1:length(neighborLabels)
                if neighborLabels(v) ~= neighbor
                    if dot( -floodGraph{label,3}(v,:), dirFromLabelToNeighbor ) < vectorsAlignedThresh
                        neighborsBorderMasked( neighborsBorderMasked == neighborLabels(v) ) = 0;
                    end
                end
            end
            for s = 1:numSamples
                center = [rows(s) cols(s)];
                localRatioPairs(s, 1) = projDistWithStencil( neighborsBorderMasked, dirFromLabelToNeighbor, center, stencil, 1 );
            end

            % compute other side of shared point with neighbor's neighbors
            currLabelPix = labeled == neighbor;
            neighbors = imdilate(currLabelPix, nHood) & ~currLabelPix;
            neighborsBorderMasked = maskImage( labeled, neighbors );
            neighborsBorderMasked( neighborsBorderMasked == label ) = 0;

            % ignore if would cause loop
            for v = 1:length(neighborNeighborLabels)
                if neighborNeighborLabels(v) ~= label
                    if dot( floodGraph{neighbor,3}(v,:), dirFromLabelToNeighbor) < vectorsAlignedThresh
                        neighborsBorderMasked( neighborsBorderMasked == neighborNeighborLabels(v) ) = 0;
                    end
                end
            end
            for s = 1:numSamples
                if localRatioPairs(s,1) ~= -1 % save some computation
                    center = [rows(s) cols(s)];
                    localRatioPairs(s, 2) = projDistWithStencil( neighborsBorderMasked, dirFromLabelToNeighbor, center, stencil, 1 );
                end
            end
            
            logRatios = zeros( numSamples, 1);
            validSample = 0;
            for s = 1:numSamples
                % if local point was able to find length on both sides
                if localRatioPairs(s,1) > 0 && localRatioPairs(s,2) > 0
                    validSample = validSample + 1;
                    logRatios(validSample) = log10( localRatioPairs(s,1) / localRatioPairs(s,2) );
                end
            end
            logRatios = logRatios(1:validSample);
            logRatio = mean(logRatios);
            distPairMap(keyPair) = logRatio;

            % check to see if candidate for Start/End label
            if abs(logRatio) > 0.5
                if logRatio > 0
                   startLabels = [startLabels; label]; 
                else
                   startLabels = [startLabels; neighbor]; 
                end
            end
        end % end IF hashmap contains pair-ratio
    end % end FOR neighbors list
end 
startLabels = unique(startLabels);

%% Compute Paths between start-end regions, Analyze Paths for pixel-lengths, ratios, and cleanup
import java.util.LinkedList;
adjMat = zeros(numlabels, numlabels);
potentialPaths = cell(1, 3);
validPathIdx = 0;
codeVote = zeros( runningVoteThreshold, codeLength );
for sEL = 1:length(startLabels)
    startLabel = startLabels(sEL);
    floodGraph{startLabel, 5} = 0;

    % convert graph connections into directed edges:
    labelToVisit = LinkedList();
    labelToVisit.add(startLabel);
    while labelToVisit.size() > 0
        label = labelToVisit.remove();
        neighborLabels = floodGraph{label, 1};
        neighboringNodes = neighborLabels;

        for n = 1:length(neighborLabels)
            neighbor = neighborLabels(n);

            % confirm direction is valid first, but ensure
            % a) all directions valid from start label, and no "prev" to
            % work with
            % b) dont revisit a path thats been analyzed already
            validDirection = 1;
            if label ~= startLabel && adjMat( neighbor, label ) == 0
                labelWasAddedByPrev = floodGraph{label, 4};
                indexOfLabelInPrev = find( floodGraph{labelWasAddedByPrev,1} == label );
                initialDir = floodGraph{labelWasAddedByPrev,3}(indexOfLabelInPrev,:);
                nextDir = floodGraph{label,3}(n,:);
                if dot( initialDir, nextDir ) < vectorsAlignedThresh
                   validDirection = 0;
                end
            end            
            
            
            if ismember(neighbor, startLabels)
               neighboringNodes( neighboringNodes == neighbor ) = []; 
            end
            
            if adjMat( neighbor, label ) == 0 && validDirection
                adjMat(label, neighbor ) = 1; %only one way directed edges allowed
                floodGraph{neighbor, 4} = label; % neighbor was added by label
                floodGraph{neighbor, 5} = floodGraph{label, 5} + 1; % number of nodes to get here
                if ~ismember(neighbor, startLabels) && floodGraph{neighbor, 5} < codeLength + 3
                    labelToVisit.add( neighbor );
                end
            else
                neighboringNodes( neighboringNodes == neighbor ) = [];
            end
        end
        floodGraph{label, 2} = neighboringNodes;
    end
    
    % traverse directed nodes for all paths from this startLabel    
    allPaths = potentialPathsSearch( floodGraph, 2, startLabel, {} );
    % cleanup paths for only those that are long enough and start/end correctly
    for potPath = 1:length(allPaths)
        if length( allPaths{potPath,1} ) == codeLength + 4
            validPathIdx = validPathIdx + 1;
            path = allPaths{potPath,1};
            potentialPaths{validPathIdx, 1} = path;

            pathCode = zeros(size(path,1) - 1,1);
            for p = 1:length(pathCode)
                
                keyPair = sprintf('%d_%d', sort([path(p,1), path(p+1,1)]));
                val = distPairMap(keyPair);
                potentialPaths{validPathIdx,2}(p,1) = val;
                if val < -0.1249 || val > 0.1761 % threshold clusters
                    potentialPaths{validPathIdx,3}(p,1) = 1;
                else
                    potentialPaths{validPathIdx,3}(p,1) = 0;
                end
            end
            
            % flip direction
            if potentialPaths{validPathIdx,3}(2,1) == 1 && potentialPaths{validPathIdx,3}(end,1) == 0
                potentialPaths{validPathIdx,1} = flipud( potentialPaths{validPathIdx,1} );
                potentialPaths{validPathIdx,2} = flipud( potentialPaths{validPathIdx,2} );
                potentialPaths{validPathIdx,3} = flipud( potentialPaths{validPathIdx,3} );
                codeVote(validPathIdx,:) = potentialPaths{validPathIdx,3}(2:end-2,1);
            else
                codeVote(validPathIdx,:) = potentialPaths{validPathIdx,3}(3:end-1,1);
            end
            
            
            if validPathIdx > runningVoteThreshold
                binCounts = zeros(codeLength, 2);
                confidence = 1.0;
                for c = 1:codeLength
                    numOnes = nnz(codeVote(1:validPathIdx,c));
                    numZeros = validPathIdx - numOnes;
                    binCounts(c,:) = [numZeros, numOnes];

                    votesForBit = max(numZeros, numOnes);
                    confidence = min( confidence, votesForBit/validPathIdx);
                end
                if confidence >= confidenceThreshold
                    electedCode = zeros(1, codeLength);
                    for c = 1:codeLength
                        if binCounts(c,1) <= binCounts(c,2)
                            electedCode(1,c) = 1;        
                        end
                    end
                    result = isequal(electedCode, correctEncoding);
%                     fprintf('elected Early (%d confidence with %d paths) Success: %d \n', confidence, validPathIdx, isequal(electedCode, correctEncoding) );
%                     fprintf(' Decoded: %s\n', sprintf('%d', electedCode))
%                     fprintf('Original: %s\n', sprintf('%d', correctEncoding ) );
%                     
                    figure; bar( binCounts );
                    xticks(1:length(binCounts));
                    text(1:length(binCounts),binCounts(:,1),num2str(binCounts(:,1)),'vert','bottom','horiz','right'); 
                    text(1:length(binCounts),binCounts(:,2),num2str(binCounts(:,2)),'vert','bottom','horiz','left');                    
                    return;
                else
                    runningVoteThreshold = runningVoteThreshold * 2;
                end
            end            
            
            
        end
    end
    

end
codeVote = codeVote(1:validPathIdx,:);

if validPathIdx == 0 && withFilter
    close all;
%     fprintf('num paths detected: none, withFilter\n');
    result = runningVote( filename, 0 );
    return;
elseif validPathIdx == 0 && ~withFilter
%     fprintf('num paths detected: none, with no Filter\n');
    clear variables;
    close all;
    result = 0;
    return;
end

fprintf('num paths detected: %d\n', size(potentialPaths,1));

%% Compare elected code with correct Encoding
binCounts = zeros(codeLength, 2);
for c = 1:codeLength
    numOnes = nnz(codeVote(:,c));
    numZeros = size(codeVote, 1) - numOnes;
    binCounts(c,:) = [numZeros, numOnes];
end
electedCode = zeros(1, codeLength);
for c = 1:codeLength
    if binCounts(c,1) <= binCounts(c,2)
        electedCode(1,c) = 1;        
    end
end

% figure; bar( binCounts );
% xticks(1:length(binCounts));
% text(1:length(binCounts),binCounts(:,1),num2str(binCounts(:,1)),'vert','bottom','horiz','right'); 
% text(1:length(binCounts),binCounts(:,2),num2str(binCounts(:,2)),'vert','bottom','horiz','left'); 

fprintf('elected Success: %d \n', isequal(electedCode, correctEncoding) );
fprintf(' Decoded: %s\n', sprintf('%d', electedCode))
fprintf('Original: %s\n', sprintf('%d', correctEncoding ) );

result = isequal(electedCode, correctEncoding);

if result
    clear variables;
    result = 1;
else
    clear variables;
    result = 0;
end
% close all;

end

