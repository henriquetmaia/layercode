function [paths] = potentialPathsSearch(floodGraph, floodGraphID, path, paths )
    regionToView = path(end);
    if ~isempty(floodGraph{regionToView,floodGraphID})
        for i = 1:length( floodGraph{regionToView,floodGraphID} )
            nextToVisit = floodGraph{regionToView,floodGraphID}(i);
            if ~ismember( nextToVisit, path )
                newPath = [path;nextToVisit];
                paths = potentialPathsSearch( floodGraph, floodGraphID, newPath, paths);               
            else
                paths{end + 1, 1} = path;
            end
        end
    else
       paths{end + 1, 1} = path;
    end
end
