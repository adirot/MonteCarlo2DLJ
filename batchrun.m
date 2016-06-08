% m = 6;
% i = 1;
% for t = [0.45 0.6 0.8 1 1.5 2] 
%     for r= [0.005 0.01 0.05 0.1] 
%         disp(i); 
%         jobs{i} = batch(c,@MC2DLJoutput,...
%             1,{625,t,r,1,'auto',10,(m/12)^(1/(m-12))/2,...
%             'pressure',true,'m',m});
%         wait(jobs{i});
%         diary(jobs{i});
%         r = fetchOutputs(jobs{i}); 
%         M(i) = r{1};
%         i = i+1;
%     end
% end
% 
% for ii = 1:length(M) 
%     jobs{ii} = batch(c,@M.MonteCarlo,1,{M(ii),10^7,1});
% end
% 
% for ii = 1:length(M) 
%     wait(jobs{ii});
%     disp(ii);
%     diary(jobs{ii});
%     r = fetchOutputs(jobs{ii}); 
%     M(ii) = r{1};
% end
% 
% for ii = 1:length(M)
%     jobs{ii} = batch(c,@M.calcRDF,1,{M(ii),10,300});
% end
% 
% for ii = 1:length(M)
%     wait(jobs{ii});
%     disp(ii);
%     diary(jobs{ii});
%     r = fetchOutputs(jobs{ii}); 
%     M(ii) = r{1};
% end

% save space - only save after 10,50,100,250,500,750,1k,2k,6k,12k steps
% (equilibration - 4000 not saved)

c = parcluster();
i = 1;
for t = [0.45 0.6 0.8 1 1.5 2] 
    for r = [0.005 0.01 0.05 0.1] 
        disp(i); 
        jobs{i} = batch(c,@MC2DLJoutput,...
            1,{625,t,r,1,'auto',10,(m/12)^(1/(m-12))/2,...
            'pressure',false,'m',m,'ufunc',@(r) 0.75*r.^-m,...
            'angleDependent',true,'angleDependence',@(a,t) 3*(cos(t).^2+cos(t-a).^2-5*cos(t).^2.*cos(t-a).^2-1/3) - (1-1)*cos(a).^2-3*(1-3)*cos(a).*cos(t).*cos(a-t),...
            'maxdAng',pi/4,'hardCoreRepRad',(m/12)^(1/(m-12))/2,...
            'fileNameInit','angdep/'});
        wait(jobs{i});
        diary(jobs{i});
        re = fetchOutputs(jobs{i}); 
        M(i) = re{1};
        i = i+1;
    end
end

for ii = 1:length(M) 
    % equilibration - 4000 steps
    temp = M(ii);
    jobs{ii} = batch(c,@temp.MonteCarlo,1,{4000,0});
    
end

for ii = 1:length(M) 
    wait(jobs{ii});
    re = fetchOutputs(jobs{ii}); 
    M(ii) = re{1};
        
    disp(['eq done in ' num2str(ii) ' of ' num2str(length(M))]);
    temp = M(ii);
    jobs{ii} = batch(c,@temp.addStep2data,1,{M(ii).currentStep,...
                        M(ii).currentCoords,...
                        M(ii).currentDists,M(ii).currentU,...
                        M(ii).currentVir,M(ii).currentPressure,...
                        M(ii).moveCount,M(ii).currentmaxdr,M(ii).Ulrc,...
                        M(ii).Plrc,...
                        M(ii).simulationParam.angleDependent,...
                        M(ii).currentAngs,...
                        M(ii).currentAlphas,M(ii).currentThetas});
    wait(jobs{ii});
    re = fetchOutputs(jobs{ii}); 
    M(ii) = re{1};
    
end


for ii = 1:length(M) 
    % run 10 more steps and save
    disp(['10 steps done in ' num2str(ii) ' of ' num2str(length(M))]);
    temp = M(ii);
    jobs{ii} = batch(c,@temp.MonteCarlo,1,{625*10,1});
    
end

for ii = 1:length(M)
    wait(jobs{ii});
    disp(ii);
    re = fetchOutputs(jobs{ii}); 
     
    M(ii) = re{1};
    
    % RDF
    newStepInd = [1 10*625];
    ind = length(newStepInd);
    
    temp = M(ii);
    jobs{ii} = batch(c,@temp.calcRDF,1,{10,300,'save2data',false,...
        'startFrom',ind-1});
    wait(jobs{ii});
    r = fetchOutputs(jobs{ii});
    RDFhisto10{ii} = mean(r{3});
    save('RDFangleDep.mat','RDFhisto10');
    % deleted unneeded info
    M(ii).data.stepInd = newStepInd;
    newCoords = M(ii).data.allCoords(:,:,1:ind);
    M(ii).data.allCoords = newCoords;
    newDists = M(ii).data.allDists(:,:,1:ind);
    M(ii).data.allDists = newDists;
    newU = M(ii).data.U(:,1:ind);
    M(ii).data.allU = newU;
    newU = M(ii).data.U(:,1:ind);
    M(ii).data.allU = newU;
    newAngs = M(ii).data.allAngs(:,:,1:ind);
    M(ii).data.allAngs = newAngs;
    newAlphas = M(ii).data.allAlphas(:,:,1:ind);
    M(ii).data.allAlphas = newAlphas;
    newThetas = M(ii).data.allThetas(:,:,1:ind);
    M(ii).data.allThetas = newThetas;
end


for ii = 1:length(M) 
    % run 100 steps and save
    
    disp(['100 steps done in ' num2str(ii) ' of ' num2str(length(M))]);
    temp = M(ii);
    jobs{ii} = batch(c,@temp.MonteCarlo,1,{625*(100-10),1});
    
end

for ii = 1:length(M)
    wait(jobs{ii});
    r = fetchOutputs(jobs{ii}); 
    M(ii) = r{1};
    
    % RDF
    newStepInd = [1 10*625 100*625];
    ind = length(newStepInd);
    
    temp = M(ii);
    jobs{ii} = batch(c,@temp.calcRDF,1,{10,300,'save2data',false,...
        'startFrom',ind-1});
    wait(jobs{ii});
    r = fetchOutputs(jobs{ii});
    RDFhistosum100{ii} = RDFhisto10{ii} + sum(r{3});
    RDFhisto100{ii} = RDFhistosum100{ii}/100;
    save('RDFangleDep.mat','RDFhisto10','RDFhisto100');
    
    % delete unneeded info
    M(ii).data.stepInd = newStepInd;
    newCoords = M(ii).data.allCoords(:,:,1:ind);
    M(ii).data.allCoords = newCoords;
    newDists = M(ii).data.allDists(:,:,1:ind);
    M(ii).data.allDists = newDists;
    newU = M(ii).data.U(:,1:ind);
    M(ii).data.allU = newU;
    newU = M(ii).data.U(:,1:ind);
    M(ii).data.allU = newU;
    newAngs = M(ii).data.allAngs(:,:,1:ind);
    M(ii).data.allAngs = newAngs;
    newAlphas = M(ii).data.allAlphas(:,:,1:ind);
    M(ii).data.allAlphas = newAlphas;
    newThetas = M(ii).data.allThetas(:,:,1:ind);
    M(ii).data.allThetas = newThetas;
end

for ii = 1:length(M) 
    % run 250 steps and save
    disp(['250 steps done in ' num2str(ii) ' of ' num2str(length(M))]);
    temp = M(ii);
    jobs{ii} = batch(c,@temp.MonteCarlo,1,{625*(250-100),1});
end

for ii = 1:length(M)
    wait(jobs{ii});
    r = fetchOutputs(jobs{ii}); 
    M(ii) = r{1};
    
    % RDF
    newStepInd = [1/625 10 100 250]*625;
    ind = length(newStepInd);
    
    temp = M(ii);
    jobs{ii} = batch(c,@temp.calcRDF,1,{10,300,'save2data',false,...
        'startFrom',ind-1});
    wait(jobs{ii});
    r = fetchOutputs(jobs{ii});
    RDFhistosum250{ii} = RDFhisto100{ii} + sum(r{3});
    RDFhisto250{ii} = RDFhistosum250{ii}/250;
    save('RDFangleDep.mat','RDFhisto10','RDFhisto100','RDFhisto250');
    
    % delete unneeded info
    M(ii).data.stepInd = newStepInd;
    newCoords = M(ii).data.allCoords(:,:,1:ind);
    M(ii).data.allCoords = newCoords;
    newDists = M(ii).data.allDists(:,:,1:ind);
    M(ii).data.allDists = newDists;
    newU = M(ii).data.U(:,1:ind);
    M(ii).data.allU = newU;
    newU = M(ii).data.U(:,1:ind);
    M(ii).data.allU = newU;
    newAngs = M(ii).data.allAngs(:,:,1:ind);
    M(ii).data.allAngs = newAngs;
    newAlphas = M(ii).data.allAlphas(:,:,1:ind);
    M(ii).data.allAlphas = newAlphas;
    newThetas = M(ii).data.allThetas(:,:,1:ind);
    M(ii).data.allThetas = newThetas;
end

for ii = 1:length(M) 
    % run 500 steps and save
    disp(['500 steps done in ' num2str(ii) ' of ' num2str(length(M))]);
    temp = M(ii);
    jobs{ii} = batch(c,@temp.MonteCarlo,1,{625*(500-250),1});
end

for ii = 1:length(M)
    wait(jobs{ii});
    r = fetchOutputs(jobs{ii}); 
    M(ii) = r{1};
    
    % RDF
    newStepInd = [1/625 10 100 250 500]*625;
    ind = length(newStepInd);
    
    temp = M(ii);
    jobs{ii} = batch(c,@temp.calcRDF,1,{10,300,'save2data',false,...
        'startFrom',ind-1});
    wait(jobs{ii});
    r = fetchOutputs(jobs{ii});
    RDFhistosum500{ii} = RDFhisto250{ii} + sum(r{3});
    RDFhisto500{ii} = RDFhistosum500{ii}/500;
    save('RDFangleDep.mat','RDFhisto10','RDFhisto100','RDFhisto250',...
        'RDFhisto500');
    
    % delete unneeded info
    M(ii).data.stepInd = newStepInd;
    newCoords = M(ii).data.allCoords(:,:,1:ind);
    M(ii).data.allCoords = newCoords;
    newDists = M(ii).data.allDists(:,:,1:ind);
    M(ii).data.allDists = newDists;
    newU = M(ii).data.U(:,1:ind);
    M(ii).data.allU = newU;
    newU = M(ii).data.U(:,1:ind);
    M(ii).data.allU = newU;
    newAngs = M(ii).data.allAngs(:,:,1:ind);
    M(ii).data.allAngs = newAngs;
    newAlphas = M(ii).data.allAlphas(:,:,1:ind);
    M(ii).data.allAlphas = newAlphas;
    newThetas = M(ii).data.allThetas(:,:,1:ind);
    M(ii).data.allThetas = newThetas;
end


for ii = 1:length(M) 
    % run 750 steps and save
    disp(['750 steps done in ' num2str(ii) ' of ' num2str(length(M))]);
    
    temp = M(ii);
    jobs{ii} = batch(c,@temp.MonteCarlo,1,{625*(750-500),1});
end

for ii = 1:length(M)
    wait(jobs{ii});
    r = fetchOutputs(jobs{ii}); 
    M(ii) = r{1};
    
    % RDF
    newStepInd = [1/625 10 100 250 500 750]*625;
    ind = length(newStepInd);
    
    temp = M(ii);
    jobs{ii} = batch(c,@temp.calcRDF,1,{10,300,'save2data',false,...
        'startFrom',ind-1});
    wait(jobs{ii});
    r = fetchOutputs(jobs{ii});
    RDFhistosum750{ii} = RDFhisto500{ii} + sum(r{3});
    RDFhisto750{ii} = RDFhistosum750{ii}/750;
    save('RDFangleDep.mat','RDFhisto10','RDFhisto100','RDFhisto250',...
        'RDFhisto500','RDFhisto750');
    
    % delete unneeded info
    M(ii).data.stepInd = newStepInd;
    newCoords = M(ii).data.allCoords(:,:,1:ind);
    M(ii).data.allCoords = newCoords;
    newDists = M(ii).data.allDists(:,:,1:ind);
    M(ii).data.allDists = newDists;
    newU = M(ii).data.U(:,1:ind);
    M(ii).data.allU = newU;
    newU = M(ii).data.U(:,1:ind);
    M(ii).data.allU = newU;
    newAngs = M(ii).data.allAngs(:,:,1:ind);
    M(ii).data.allAngs = newAngs;
    newAlphas = M(ii).data.allAlphas(:,:,1:ind);
    M(ii).data.allAlphas = newAlphas;
    newThetas = M(ii).data.allThetas(:,:,1:ind);
    M(ii).data.allThetas = newThetas;
end


for ii = 1:length(M) 
    % run 1000 steps and save
    disp(['1000 steps done in ' num2str(ii) ' of ' num2str(length(M))]);
    temp = M(ii);
    jobs{ii} = batch(c,@temp.MonteCarlo,1,{625*(1000-750),1});
end

for ii = 1:length(M)
    wait(jobs{ii});
    r = fetchOutputs(jobs{ii}); 
    M(ii) = r{1};
    
    % RDF
    newStepInd = [1/625 10 100 250 500 750 1000]*625;
    ind = length(newStepInd);
    
    temp = M(ii);
    jobs{ii} = batch(c,@temp.calcRDF,1,{10,300,'save2data',false,...
        'startFrom',ind-1});
    wait(jobs{ii});
    r = fetchOutputs(jobs{ii});
    RDFhistosum1000{ii} = RDFhisto750{ii} + sum(r{3});
    RDFhisto1000{ii} = RDFhistosum1000{ii}/1000;
    save('RDFangleDep.mat','RDFhisto10','RDFhisto100','RDFhisto250',...
        'RDFhisto500','RDFhisto750','RDFhisto1000');
    
    % delete unneeded info
    M(ii).data.stepInd = newStepInd;
    newCoords = M(ii).data.allCoords(:,:,1:ind);
    M(ii).data.allCoords = newCoords;
    newDists = M(ii).data.allDists(:,:,1:ind);
    M(ii).data.allDists = newDists;
    newU = M(ii).data.U(:,1:ind);
    M(ii).data.allU = newU;
    newU = M(ii).data.U(:,1:ind);
    M(ii).data.allU = newU;
    newAngs = M(ii).data.allAngs(:,:,1:ind);
    M(ii).data.allAngs = newAngs;
    newAlphas = M(ii).data.allAlphas(:,:,1:ind);
    M(ii).data.allAlphas = newAlphas;
    newThetas = M(ii).data.allThetas(:,:,1:ind);
    M(ii).data.allThetas = newThetas;
end


for ii = 1:length(M) 
    % run 2000 steps and save
    disp(['2000 steps done in ' num2str(ii) ' of ' num2str(length(M))]);
    temp = M(ii);
    jobs{ii} = batch(c,@temp.MonteCarlo,1,{625*(2000-1000),1});
end

for ii = 1:length(M)

    wait(jobs{ii});
    r = fetchOutputs(jobs{ii}); 
    M(ii) = r{1};
    
    % RDF
    newStepInd = [1/625 10 100 250 500 750 1000 2000]*625;
    ind = length(newStepInd);
    
    temp = M(ii);
    jobs{ii} = batch(c,@temp.calcRDF,1,{10,300,'save2data',false,...
        'startFrom',ind-1});
    wait(jobs{ii});
    r = fetchOutputs(jobs{ii});
    RDFhistosum2000{ii} = RDFhisto1000{ii} + sum(r{3});
    RDFhisto2000{ii} = RDFhistosum2000{ii}/2000;
    save('RDFangleDep.mat','RDFhisto10','RDFhisto100','RDFhisto250',...
        'RDFhisto500','RDFhisto750','RDFhisto1000','RDFhisto2000');
    
    % delete unneeded info
    M(ii).data.stepInd = newStepInd;
    newCoords = M(ii).data.allCoords(:,:,1:ind);
    M(ii).data.allCoords = newCoords;
    newDists = M(ii).data.allDists(:,:,1:ind);
    M(ii).data.allDists = newDists;
    newU = M(ii).data.U(:,1:ind);
    M(ii).data.allU = newU;
    newU = M(ii).data.U(:,1:ind);
    M(ii).data.allU = newU;
    newAngs = M(ii).data.allAngs(:,:,1:ind);
    M(ii).data.allAngs = newAngs;
    newAlphas = M(ii).data.allAlphas(:,:,1:ind);
    M(ii).data.allAlphas = newAlphas;
    newThetas = M(ii).data.allThetas(:,:,1:ind);
    M(ii).data.allThetas = newThetas;
end


for ii = 1:length(M) 
    % run 6000 steps and save
    disp(['6000 steps done in ' num2str(ii) ' of ' num2str(length(M))]);
    
    temp = M(ii);
    jobs{ii} = batch(c,@temp.MonteCarlo,1,{625*(6000-2000),1});
    wait(jobs{ii});
    
end

for ii = 1:length(M)

    r = fetchOutputs(jobs{ii}); 
    M(ii) = r{1};
    
    % RDF
    newStepInd = [1/625 10 100 250 500 750 1000 2000 6000]*625;
    ind = length(newStepInd);
    
    temp = M(ii);
    jobs{ii} = batch(c,@temp.calcRDF,1,{10,300,'save2data',false,...
        'startFrom',ind-1});
    wait(jobs{ii});
    r = fetchOutputs(jobs{ii});
    RDFhistosum6000{ii} = RDFhisto2000{ii} + sum(r{3});
    RDFhisto6000{ii} = RDFhistosum6000{ii}/6000;
    save('RDFangleDep.mat','RDFhisto10','RDFhisto100','RDFhisto250',...
        'RDFhisto500','RDFhisto750','RDFhisto1000','RDFhisto2000',...
        'RDFhisto6000');
    
    % delete unneeded info
    M(ii).data.stepInd = newStepInd;
    newCoords = M(ii).data.allCoords(:,:,1:ind);
    M(ii).data.allCoords = newCoords;
    newDists = M(ii).data.allDists(:,:,1:ind);
    M(ii).data.allDists = newDists;
    newU = M(ii).data.U(:,1:ind);
    M(ii).data.allU = newU;
    newU = M(ii).data.U(:,1:ind);
    M(ii).data.allU = newU;
    newAngs = M(ii).data.allAngs(:,:,1:ind);
    M(ii).data.allAngs = newAngs;
    newAlphas = M(ii).data.allAlphas(:,:,1:ind);
    M(ii).data.allAlphas = newAlphas;
    newThetas = M(ii).data.allThetas(:,:,1:ind);
    M(ii).data.allThetas = newThetas;
end


for ii = 1:length(M) 
    % run 12k steps and save
    disp(['12k steps done in ' num2str(ii) ' of ' num2str(length(M))]);
    temp = M(ii);
    jobs{ii} = batch(c,@temp.MonteCarlo,1,{625*(12000-6000),1});
    
end

for ii = 1:length(M)

    wait(jobs{ii});
    r = fetchOutputs(jobs{ii}); 
    M(ii) = r{1};
    
    % RDF
    newStepInd = [1/625 10 100 250 500 750 1000 2000 6000 12000]*625;
    ind = length(newStepInd);
    
    temp = M(ii);
    jobs{ii} = batch(c,@temp.calcRDF,1,{10,300,'save2data',false,...
        'startFrom',ind-1});
    wait(jobs{ii});
    r = fetchOutputs(jobs{ii});
    RDFhistosum12000{ii} = RDFhisto6000{ii} + sum(r{3});
    RDFhisto12000{ii} = RDFhistosum12000{ii}/12000;
    save('RDFangleDep.mat','RDFhisto10','RDFhisto100','RDFhisto250',...
        'RDFhisto500','RDFhisto750','RDFhisto1000','RDFhisto2000',...
        'RDFhisto6000','RDFhisto12000');
    
    % delete unneeded info
    M(ii).data.stepInd = newStepInd;
    newCoords = M(ii).data.allCoords(:,:,1:ind);
    M(ii).data.allCoords = newCoords;
    newDists = M(ii).data.allDists(:,:,1:ind);
    M(ii).data.allDists = newDists;
    newU = M(ii).data.U(:,1:ind);
    M(ii).data.allU = newU;
    newU = M(ii).data.U(:,1:ind);
    M(ii).data.allU = newU;
    newAngs = M(ii).data.allAngs(:,:,1:ind);
    M(ii).data.allAngs = newAngs;
    newAlphas = M(ii).data.allAlphas(:,:,1:ind);
    M(ii).data.allAlphas = newAlphas;
    newThetas = M(ii).data.allThetas(:,:,1:ind);
    M(ii).data.allThetas = newThetas;
end

