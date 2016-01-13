function ind = canGetUfromgR(isotherms,fit,maxres)
    
    for i = length(isotherms)
        Plin = fit{1,1}(1,i).P;
        Pvir = fit{1,2}(1,i).P;
        Pdata = isotherms(i).pressure;
        
        reslin = abs((Plin-Pdata)./Pdata);
        resvir = abs((Pvir-Pdata)./Pdata);
        
        ind{1,i} = find(and((reslin >= maxres),(resvir <= maxres)));
        
    end
end
