for i = 1:length(Outm1) 
    for j = 1:length(leg) 
        if strcmp(leg{1,j},Outm1{1,i}) 
            keep1(i) = 1; 
            break;
        else
            keep1(i) = 0;
        end; 
    end; 
end;