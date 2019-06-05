function projDist = projDistWithStencil( otherMap, projDir, startingPoint, stencil, useMin )

    topLeft = floor(startingPoint - size(stencil)/2 );
    imagePatch = otherMap(topLeft(1):topLeft(1)+size(stencil,1)-1,topLeft(2):topLeft(2)+size(stencil,2)-1);
    imagePatch = single( imagePatch ~= 0 ); % remove label bias by binarization, but then make single for mult
    distToStartingPoint = times(imagePatch, stencil);
    tmp = (distToStartingPoint == 0);
    midPt = ceil(size(stencil)/2);
    if useMin
        distToStartingPoint(tmp) = Inf;
        [~, indMin] = min(distToStartingPoint(:));
        [nr,nc] = ind2sub( size(distToStartingPoint), indMin);
%         [nr,nc] = find( distToStartingPoint == min(distToStartingPoint(:)));
%         nr = nr(1);
%         nc = nc(1);
        if distToStartingPoint(nr,nc) == Inf
            projDist = -1;
            return;
        end
        projDist = abs(dot(( [nr,nc] - midPt ), projDir));
        
    else % if using max, should only form line on near side...
        distToStartingPoint(tmp) = -Inf;
        [nr,nc] = find( distToStartingPoint == max(distToStartingPoint(:)));
        nr = nr(1);
        nc = nc(1);        
        if distToStartingPoint(nr,nc) == -Inf
            projDist = -1;
            return;
        end
        projDist = abs(dot(( [nr,nc] - midPt ), projDir));
    end
end
