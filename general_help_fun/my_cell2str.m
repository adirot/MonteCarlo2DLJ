function string = my_cell2str(cellWithstrs)
    string = '';
    for i = 1:length(cellWithstrs)
        string = [string,' ',cellWithstrs{i}];
    end
    string = string(2:end);
    
end