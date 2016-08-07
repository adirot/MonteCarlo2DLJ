function newLogFileName = createLogFile(logFileInit,oldLogFileName,totSec)
    
if ~isempty(oldLogFileName)
    delete(oldLogFileName);
end

newLogFileName = [logFileInit 'secPassed'...
    my_num2str(totSec) 'date' nowdatetimestr() '.txt'];
fileID = fopen(newLogFileName, 'w');
fclose(fileID);

end