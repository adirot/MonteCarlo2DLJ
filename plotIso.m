function plotIso(isotherms)
    
    Niso = length(isotherms);
    Nrho = length(isotherms(1).rho);
    pressure = zeros(Nrho,Niso);
    rho = zeros(Nrho,Niso);
    
    for i = 1:Niso
        rho(:,i) = isotherms(i).rho;
        pressure(:,i) = isotherms(i).pressure;
        leg{1,i} = ['T = ' num2str(isotherms(i).T)];
        
    end
    
    colorPlot(rho,pressure,'addLegend',leg);
end