m = 6;i = 1; for t = [0.45 0.6 0.8 1 1.5 2] for r= [0.005 0.01 0.05 0.1] disp(i); jobs{i} = batch(c,@MC2DLJoutput,1,{625,t,r,1,'auto',10,(m/12)^(1/(m-12))/2,'pressure',true,'m',m});wait(jobs{i});diary(jobs{i});r = fetchOutputs(jobs{i}); M(i) = r{1};i = i+1;end;end; 
for ii = 1:length(M) jobs{ii} = batch(c,@MonteCarlo,1,{M(ii),10^7,1});end;
for ii = 1:length(M) wait(jobs{ii});disp(ii);diary(jobs{ii});r = fetchOutputs(jobs{ii}); M(ii) = r{1};end;
for ii = 1:length(M); jobs{ii} = batch(c,@calcRDF,1,{M(ii),10,300});end; 
for ii = 1:length(M) wait(jobs{ii});disp(ii);diary(jobs{ii});r = fetchOutputs(jobs{ii}); M(ii) = r{1};end;
