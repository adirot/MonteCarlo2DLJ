function exists = existInMatfile(fileName,varName)
    var = whos('-file',fileName);
    var = {var.name};
    exists = any(strcmp(var,varName));
end