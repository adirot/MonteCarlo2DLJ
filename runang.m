function M = runang(t,r,m)
v = 1/2; % poisson ratio 
a2 = (1-v)/v;
M = MC2DLJoutput(2,t,r,1,'auto',10,(m/12)^(1/(m-12))/2,...
            'pressure',false,'m',m,'ufunc',@(r) r.^-m,...
            'angleDependent',true,...
            'angleDependence',@(a,t) 3*(cos(t).^2+cos(t-a).^2-5*cos(t).^2.*cos(t-a).^2-1/3) - (1-a2)*cos(a).^2-3*(a2-3)*cos(a).*cos(t).*cos(a-t),...
            'maxdAng',pi/4,'hcr',true,...
            'fileNameInit','angdep/','dontSaveDists',true);
        
% equilibration
tic;
M = M.MonteCarlo(8000,1);
toc
disp(['done equilibrating ' num2str(t) ' ' num2str(r) ' ' num2str(m)]);
M = M.addStep2data(M.currentSweep,...
                        M.currentCoords,...
                        M.currentDists,M.currentU,...
                        M.currentVir,M.currentPressure,...
                        M.moveCount,M.currentmaxdr,M.Ulrc,...
                        M.Plrc,...
                        M.simulationParam.angleDependent,...
                        M.currentAngs,...
                        M.currentAlphas,M.currentThetas,...
                        M.simulationParam.dontSaveDists);
      
% run
disp(['starting run - ' num2str(t) ' ' num2str(r) ' ' num2str(m)]);
tic;

%% run on 10 differet workers:
for i = 1:10
    Mcopies(i) = M.copyOutputFile(num2str(i));
end

% separate MC output file to 10 different ones:

c = parcluster();
for i = 1:10
        temp = Mcopies(i);
        jobs{i} = batch(c,@temp.MonteCarlo,1,{16000/10,10});
end

for i = 1:10
    wait(jobs{i});
    diary(jobs{i});
    re = fetchOutputs(jobs{i}); 
    Mcopies(i) = re{1};
end

% unify runs
M = unifyMCobj(Mcopies);

% delete Mcopies
for i = 1:10
    delete(Mcopies(i).fileName);
end


toc
disp(['160000 sweeps done - ' num2str(t) ' ' num2str(r) ' ' num2str(m)]);
                    
end