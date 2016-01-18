function ind = canGetUfromgR(isotherms,fit,maxres)
    
    ind = zeros(length(isotherms),length(isotherms(1).rho));

    for i = 1:length(isotherms)
        Plin = fit{1,1}(1,i).P;
        Pvir = fit{1,2}(1,i).P;
        Pdata = isotherms(i).pressure;
        
        reslin = abs((Plin-Pdata)./Pdata);
        resvir = abs((Pvir-Pdata)./Pdata);
        
        ind(i,:) = and((reslin >= maxres),(resvir <= maxres));
        
    end
end
