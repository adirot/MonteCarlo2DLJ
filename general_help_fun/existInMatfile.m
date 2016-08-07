function exists = existInMatfile(fileNameOrMatfile,varName)
    if isobject(fileNameOrMatfile)
        try
            fileName = fileNameOrMatfile.Properties.Source;
        catch
            error('first input must be file name or matfile object');
        end
    else
        if ischar(fileNameOrMatfile)
            fileName = fileNameOrMatfile;
        else
            error('first input must be file name or matfile object');
        end
    end
    
    var = whos('-file',fileName);
    var = {var.name};
    exists = any(strcmp(var,varName));
end