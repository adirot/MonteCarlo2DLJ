function MonteCarlo2DLJ(N,Nsteps,varargin)

%% Monte-Carlo in NVT ensemble for Lennard-Jonse potantioal in 2D %%

p = inputParser();
addOptional(p, 'T', [0.01 0.1 0.2 0.3 0.45 0.6 0.7 1 4.5 10 100]);
addOptional(p, 'm', 6);
parse(p, varargin{:});
Results = p.Results;
T = Results.T;
m = Results.m;
firstSteps2ignore = 10;
        
%rhorand = [0.0001 0.0003 0.0005 0.0007 0.001 0.003 0.005 0.007 0.01 0.03 0.05 0.07 0.1 0.2 0.3 0.4];
%rhohex  = 0.5:0.1:0.7;
%rho = [rhorand, rhohex];
rho = 0.1:0.1:0.7;
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
list = {};
for i = 1:length(T)
    
    if T(i) < 10
        rcutoff = 6;
    end
    
    if T(i) == 10
        rcutoff = 4;
    end 
    
    if T(i) == 100
        rcutoff = 3;
    end 
    
    I(i) = isotherm('N',N,'T',T(i),'rho',rho,'initialmaxdr',maxdr,...
        'initialConfig','auto','rCutoff',rcutoff,'r',r,...
        'cutEquilirization',false,'m',m);
    list(i,:) = I(i).datafileList;
end



disp(['created Isotherms N = ' num2str(N) ' m = ' num2str(m)]);
save(['isoObjN' num2str(N) 'm' num2str(m) 'sigma3_6'],'-v7.3');

for i = 1:length(I)
    I(i) = I(i).calcIso(Nsteps,10);
end

disp(['calculated Isotherms N = ' num2str(N) ' m = ' num2str(m)]);
save(['isoObjN' num2str(N) 'm' num2str(m) 'sigma3_6'],'-v7.3');

for i = 1:length(I)
    I(i) = I(i).calcMeanWithoutFirstSteps(50);
    I(i) = I(i).calcCv();
end

disp(['calculated Cv N = ' num2str(N) ' m = ' num2str(m)]);
save(['isoObjN' num2str(N) 'm' num2str(m) 'sigma3_6'],'-v7.3');

for i = 1:length(T)
    
    [I(i), h] = I(i).plotPropVsStep('U');
    saveas(h,['UvsStepN' num2str(N) 'T' my_num2str(T(i)) 'm' num2str(m) 'sigma3_6' '.fig']);
    saveas(h,['UvsStepN' num2str(N) 'T' my_num2str(T(i)) 'm' num2str(m) 'sigma3_6' '.jpg']);
    close all;
    
    [I(i), h] = I(i).plotPropVsStep('P');
    saveas(h,['PvsStepN' num2str(N) 'T' my_num2str(T(i)) 'm' num2str(m) 'sigma3_6' '.fig']);
    saveas(h,['PvsStepN' num2str(N) 'T' my_num2str(T(i)) 'm' num2str(m) 'sigma3_6' '.jpg']);
    close all;

end

disp(['ploted PvsStep N = ' num2str(N) ' m = ' num2str(m)]);

[isotherms,fit,canGetUfromgRind,~,~,P,U,T,Z,Zx] = plotIso('isotherms',I,...
    'N',N,'fitprop',{'plotLin', 'plotVirialExp'},...
    'residuals',{false, false},'talk',true,'fileNameEnd','sigma3_6');

disp(['ploted isotherms N = ' num2str(N) ' m = ' num2str(m)]);
save(['isoObjN' num2str(N) 'm' num2str(m) 'sigma3_6'],'-v7.3');


clear I fit canGetUfromgRind;


%list = dir(['N' num2str(N) 'T*m' num2str(m) '*mat']);
%list = {list.name};

RDFm5 = RDFoutput('dataFileList',list,'N',N);
save(['RDFObjN' num2str(N) 'm' num2str(m) 'sigma3_6'],'-v7.3');

disp(['created RDFobj N = ' num2str(N) ' m = ' num2str(m)]);

RDFm5 = RDFm5.calcAllRDF(10,300,'skipExisting',true);
save(['RDFObjN' num2str(N) 'm' num2str(m) 'sigma3_6'],'-v7.3');

disp(['calculated RDFs N = ' num2str(N) ' m = ' num2str(m)]);

RDFm5 = RDFm5.plotRDFT('keepFigOpen',false,'fileNameEnd','sigma3_6');
save(['RDFObjN' num2str(N) 'm' num2str(m) 'sigma3_6'],'-v7.3');

RDFm5 = RDFm5.plotRDFT('keepFigOpen',false,'plotLog',true,'fileNameEnd','sigma3_6');
save(['RDFObjN' num2str(N) 'm' num2str(m) 'sigma3_6'],'-v7.3');

RDFm5 = RDFm5.plotRDFrho('keepFigOpen',false,'fileNameEnd','sigma3_6');
save(['RDFObjN' num2str(N) 'm' num2str(m) 'sigma3_6'],'-v7.3');

RDFm5 = RDFm5.plotRDFrho('keepFigOpen',false,'plotLog',true,'fileNameEnd','sigma3_6');
save(['RDFObjN' num2str(N) 'm' num2str(m) 'sigma3_6'],'-v7.3');

disp(['ploted RDF N = ' num2str(N) ' m = ' num2str(m)]);

clear RDF;

% sq = (4:10).^2;
% 
% for i = 1:length(sq)
%     RhoDistrib(i) = RhoDistriboutput('N',N);
% end
% 
% disp(['created rhoDistrib Obj N = ' num2str(N)]);
% save(['rhoDistObjN' num2str(N)],'-v7.3');
% 
% for i = 1:length(sq)
%     RhoDistrib(i) = RhoDistrib(i).calcAllrhoDistrib(sq(i),30);
%     RhoDistrib(i) = RhoDistrib(i).plotrhoDistribT('addStr2title',...
%         ['squares = ' num2str(sq(i))],...
%         'addFileNameEnd',['Sq' num2str(sq(i))],'keepFigOpen',false);
%     RhoDistrib(i) = RhoDistrib(i).plotrhoDistribrho('addStr2title',...
%         ['squares = ' num2str(sq(i))],...
%         'addFileNameEnd',['Sq' num2str(sq(i))],'keepFigOpen',false);
% end
% 
% disp(['ploted rhoDistrib Obj N = ' num2str(N)]);

end




