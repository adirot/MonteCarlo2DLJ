function MonteCarlo2DLJ(N,Nsteps)

%% Monte-Carlo in NVT ensemble for Lennard-Jonse potantioal in 2D %%


        T = [0.1 0.2 0.3 0.4 0.45 0.5 0.6 0.7 0.9 1 10 100];
        rhorand = [0.0001 0.0003 0.0005 0.0007 0.001 0.003 0.005 0.007 0.01 0.03 0.05 0.07 0.1 0.2 0.3 0.4];
        rhohex  = 0.42:0.02:0.7;
        maxdr = 1;
        rCutoff = 2.5;
        saveEvery = 10;
%         T = [1 2 3];
%         rhorand = [0.1 0.2 0.3 0.4];
%         rhohex  = [0.5 0.6 0.7];
%         maxdr = 1;
%         rCutoff = 2.5;
%         saveEvery = 1;


r = 2^(1/6)/2; % particle radius in reduced units 

for i = 1:length(T)
    Irand(i) = isotherm('N',N,'T',T(i),'rho',rhorand,'initialmaxdr',maxdr,...
        'initialConfig','random','rCutoff',2.5,'r',r,...
        'cutEquilirization',false);
    Ihex(i) = isotherm('N',N,'T',T(i),'rho',rhohex,'initialmaxdr',maxdr,...
        'initialConfig','hex','rCutoff',2.5,'r',r,...
        'cutEquilirization',false);
end

save(['isoN' num2str(N)],'-v7.3');

for i = 1:length(T)
    Irand(i) = Irand(i).calcIso(Nsteps,10);
    Ihex(i) = Ihex(i).calcIso(Nsteps,10);
end

save(['isoObjN' num2str(N)],'-v7.3');

for i = 1:length(T)
    
    [Irand(i), h] = Irand(i).plotPropVsStep('U');
    [Ihex(i), h] = Ihex(i).plotPropVsStep('U','figHandle',h);
    saveas(h,['UvsStepN' num2str(N) 'T' my_num2str(T(i)) '.fig']);
    saveas(h,['UvsStepN' num2str(N) 'T' my_num2str(T(i)) '.jpg']);
    close all;
    
    [Irand(i), h] = Irand(i).plotPropVsStep('P');
    [Ihex(i), h] = Ihex(i).plotPropVsStep('P','figHandle',h);
    saveas(h,['PvsStepN' num2str(N) 'T' my_num2str(T(i)) '.fig']);
    saveas(h,['PvsStepN' num2str(N) 'T' my_num2str(T(i)) '.jpg']);
    close all;

end

[~,fitrand,canGetUfromgRindrand,hPvsRho,hPvsV] = plotIso('isotherms',Irand,...
    'N',N,'fitprop',{'plotLin', 'plotVirialExp'},...
    'residuals',{false, false});

[~,fithex,canGetUfromgRindhex,~,~] = plotIso('isotherms',Ihex,...
    'N',N,'fitprop',{'plotLin', 'plotVirialExp'}...
    ,'hPvsRho',hPvsRho,'hPvsV',hPvsV,...
    'residuals',{false, false});

save(['isoObjN' num2str(N)],'-v7.3');

clear Irand Ihex fithex fitrand canGetUfromgRindrand canGetUfromgRindhex;

RDF = RDFoutput('N',N);
save(['RDFObjN' num2str(N)],'-v7.3');

RDF = RDF.calcAllRDF(10,300);
save(['RDFObjN' num2str(N)],'-v7.3');

RDF = RDF.plotRDFT('keepFigOpen',false);
save(['RDFObjN' num2str(N)],'-v7.3');

RDF = RDF.plotRDFT('keepFigOpen',false,'plotLog',true);
save(['RDFObjN' num2str(N)],'-v7.3');

RDF = RDF.plotRDFrho('keepFigOpen',false);
save(['RDFObjN' num2str(N)],'-v7.3');

RDF = RDF.plotRDFrho('keepFigOpen',false,'plotLog',true);
save(['RDFObjN' num2str(N)],'-v7.3');

clear RDF;

sq = (4:10).^2;

for i = 1:length(sq)
    RhoDistrib(i) = RhoDistriboutput('N',N);
end

save(['rhoDistObjN' num2str(N)],'-v7.3');

for i = 1:length(sq)
    RhoDistrib(i) = RhoDistrib(i).calcAllrhoDistrib(sq(i),30);
    RhoDistrib(i) = RhoDistrib(i).plotrhoDistribT('addStr2title',...
        ['squares = ' num2str(sq(i))],...
        'addFileNameEnd',['Sq' num2str(sq(i))],'keepFigOpen',false);
    RhoDistrib(i) = RhoDistrib(i).plotrhoDistribrho('addStr2title',...
        ['squares = ' num2str(sq(i))],...
        'addFileNameEnd',['Sq' num2str(sq(i))],'keepFigOpen',false);
end


end




