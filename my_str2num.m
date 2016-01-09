function num = my_str2num(str)
    str(str == '_') = '.';
    num = str2num(str);
end