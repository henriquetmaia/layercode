


function [] = viewAngleCheck( obj_basename )

    logFileName = sprintf('%s/statLog.txt', obj_basename );
    fileID = fopen(logFileName,'a');
    picturesName = sprintf('%s/*.png', obj_basename );
    files = dir(picturesName);
    finishedFolderName = sprintf('%s/finished_%s', obj_basename, obj_basename );
    if ~exist(finishedFolderName, 'dir')
        % Folder does not exist so create it.
        mkdir( finishedFolderName );
    end    
    
    for file = files'
        currPic = sprintf('%s/%s', obj_basename, file.name);
        if runningVote( currPic, 1 )
            fprintf(fileID,'%s %d\n',file.name, 1 );
        else
            fprintf(fileID,'%s %d\n',file.name, 0 );
        end
        movefile( currPic, finishedFolderName );
    end
    fclose(fileID);
%     exit(0);
end
