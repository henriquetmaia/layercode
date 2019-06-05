
% searches for closest pixel neighbor 
% from 'startingPoint'
% to 'colorToFind'
% using the given 'stencil' as a length measurement and
function nearestNeighbor = findWithStencil( origMap, startingPoint, colorToFind, stencil )
    nearestNeighbor = [-1, -1];
    
    topLeft = floor(startingPoint - size(stencil)/2 );
    imagePatch = origMap(topLeft(1):topLeft(1)+size(stencil,1)-1,topLeft(2):topLeft(2)+size(stencil,2)-1);

    masked = imagePatch == colorToFind;
    if max(masked(:)) < 1
        return
    end
    masked = single(masked);
    distToStartingPoint = times(masked, stencil);
    tmp = (distToStartingPoint == 0);
    distToStartingPoint(tmp) = Inf;
    
    [nr,nc] = find( distToStartingPoint == min(distToStartingPoint(:)));
    nearestNeighbor = [nr(1),nc(1)] + topLeft;
   
end
