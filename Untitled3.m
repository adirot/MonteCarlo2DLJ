
% fit all matfiles to the potantail (logRDF) for different steps

steps = {'10';'50'; '100'; '1k'; '2k'; '6k'; '12k'};
list = dir('*mat');
list = {list.name}';
start = 4000;

for i = 1:length(list)
        
        
        M = MC2DLJoutput(list{i,1});
        T = M.simulationParam.T;
        rho = M.simulationParam.rho;
        m = M.simulationParam.m;
        
        if i == 1
            x = M.data.RDFbins;
            xs = [x;x;x;x;x;x;x];
        end
        
        RDF10 = mean(M.data.RDFhisto(start:(start+10),1:300));
        RDF50 = mean(M.data.RDFhisto(start:(start+50),1:300));
        RDF100 = mean(M.data.RDFhisto(start:(start+100),1:300));
        RDF1k = mean(M.data.RDFhisto(start:(start+1000),1:300));
        RDF2k = mean(M.data.RDFhisto(start:(start+2000),1:300));
        RDF6k = mean(M.data.RDFhisto(start:(start+6000),1:300));
        RDF12k = mean(M.data.RDFhisto(start:(start+12000),1:300));
        
        ys = [-log(RDF10);-log(RDF50);-log(RDF100);-log(RDF1k);...
            -log(RDF2k);-log(RDF6k);-log(RDF12k)];
        [fitresult, gof] = createFitRDF(xs, ys, T, rho,steps);
        saveas(gcf,['fitRDFT' my_num2str(T)...
            'rho' my_num2str(rho) 'm' num2str(m) 'steps.fig']);
        saveas(gcf,['fitRDFT' my_num2str(T)...
            'rho' my_num2str(rho) 'm' num2str(m) 'steps.jpg']);
        close all;
        disp(i);
end