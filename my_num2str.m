function str = my_num2str(num)
 % the number will be converted into a string.
 % if there is a dot in the number, it will be turned into _
    
    str = num2str(num);
    dotPlace = find('.' == str);
    if ~isempty(dotPlace)
        str(dotPlace) = '_';
    end
    
end