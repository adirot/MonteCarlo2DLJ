%% get min / max values of RDF %%


steps = {'10';'50'; '100'; '1k'; '2k'; '6k'; '12k'};
list = dir('N*mat');
list = {list.name}';
start = 4000;

tind = 0;
rind = 0;
mind = 0;

for t = [0.45,0.6,0.8,1,1.5,2]
    tind = tind + 1;
    rind = 0;
    for r = [0.005,0.01,0.05,0.1]
        rind = rind + 1;
        mind = 0;
        for m = [3,4,5,6]
            mind = mind + 1; 
            
            disp(t); 
            disp(r);
            disp(m);
            disp('-----');
            
            list = dir(['*T' my_num2str(t) 'rho' my_num2str(r) '*_m' num2str(m)...
                '*mat']);
            list = {list.name}';

            if length(list) > 1
                disp('more than one file!');
            end

            M = MC2DLJoutput(list{1,1});

            if and(and(t == 0.45,r == 0.005), m==3)
                x = M.data.RDFbins;
            end
            tic;     
        
            RDF10min{tind,rind,mind} = min(M.data.RDFhisto(start:(start+10),1:300));
            RDF50min{tind,rind,mind} = min(M.data.RDFhisto(start:(start+50),1:300));
            RDF100min{tind,rind,mind} = min(M.data.RDFhisto(start:(start+100),1:300));
            RDF1kmin{tind,rind,mind} = min(M.data.RDFhisto(start:(start+1000),1:300));
            RDF2kmin{tind,rind,mind} = min(M.data.RDFhisto(start:(start+2000),1:300));
            RDF6kmin{tind,rind,mind} = min(M.data.RDFhisto(start:(start+6000),1:300));
            RDF12kmin{tind,rind,mind} = min(M.data.RDFhisto(start:(start+12000),1:300));
        
            
            RDF10max{tind,rind,mind} = max(M.data.RDFhisto(start:(start+10),1:300));
            RDF50max{tind,rind,mind} = max(M.data.RDFhisto(start:(start+50),1:300));
            RDF100max{tind,rind,mind} = max(M.data.RDFhisto(start:(start+100),1:300));
            RDF1kmax{tind,rind,mind} = max(M.data.RDFhisto(start:(start+1000),1:300));
            RDF2kmax{tind,rind,mind} = max(M.data.RDFhisto(start:(start+2000),1:300));
            RDF6kmax{tind,rind,mind} = max(M.data.RDFhisto(start:(start+6000),1:300));
            RDF12kmax{tind,rind,mind} = max(M.data.RDFhisto(start:(start+12000),1:300));
            
       
%             plot(x,RDF12k{tind,rind,mind});            
%             hold on;
%             plot(x,RDF12k{tind,rind,mind}+...
%                 (RDF12k{tind,rind,mind}.*RDF12kstd{tind,rind,mind})/sqrt(12000),'-r');
%             
%             plot(x,RDF12k{tind,rind,mind}-...
%                 (RDF12k{tind,rind,mind}.*RDF12kstd{tind,rind,mind})/sqrt(12000),'-r');
%             
%             plot(x,RDF12kmin{tind,rind,mind},'r');
%             plot(x,RDF12kmax{tind,rind,mind},'r');
%             
%             saveas(gcf,['RDFstdminmax12kT' my_num2str(t)...
%                 'rho' my_num2str(r) 'm' num2str(m) '.fig']);
%             saveas(gcf,['RDFstdminmax12kT' my_num2str(t)...
%                 'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
%             close all;
% 
%             
%             plot(x,RDF10{tind,rind,mind});            
%             hold on;
%             plot(x,RDF10{tind,rind,mind}+...
%                 (RDF10{tind,rind,mind}.*RDF10std{tind,rind,mind})/sqrt(10),'-r');
%             
%             plot(x,RDF10{tind,rind,mind}-...
%                 (RDF10{tind,rind,mind}.*RDF10std{tind,rind,mind})/sqrt(10),'-r');
%             
%             plot(x,RDF10min{tind,rind,mind},'r');
%             plot(x,RDF10max{tind,rind,mind},'r');
%             
%             saveas(gcf,['RDFstdminmax10T' my_num2str(t)...
%                 'rho' my_num2str(r) 'm' num2str(m) '.fig']);
%             saveas(gcf,['RDFstdminmax10T' my_num2str(t)...
%                 'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
%             close all;
%             
%             
%             plot(x,RDF50{tind,rind,mind});            
%             hold on;
%             plot(x,RDF50{tind,rind,mind}+...
%                 (RDF50{tind,rind,mind}.*RDF50std{tind,rind,mind})/sqrt(50),'-r');
%             
%             plot(x,RDF50{tind,rind,mind}-...
%                 (RDF50{tind,rind,mind}.*RDF50std{tind,rind,mind})/sqrt(50),'-r');
%             
%             
%             plot(x,RDF50min{tind,rind,mind},'r');
%             plot(x,RDF50max{tind,rind,mind},'r');
%             
%             saveas(gcf,['RDFstdminmax50T' my_num2str(t)...
%                 'rho' my_num2str(r) 'm' num2str(m) '.fig']);
%             saveas(gcf,['RDFstdminmax50T' my_num2str(t)...
%                 'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
%             close all;
%             
%             
%             plot(x,RDF100{tind,rind,mind});            
%             hold on;
%             plot(x,RDF100{tind,rind,mind}+...
%                 (RDF100{tind,rind,mind}.*RDF100std{tind,rind,mind})/sqrt(100),'-r');
%             
%             plot(x,RDF100{tind,rind,mind}-...
%                 (RDF100{tind,rind,mind}.*RDF100std{tind,rind,mind})/sqrt(100),'-r');
%             
%             
%             plot(x,RDF100min{tind,rind,mind},'r');
%             plot(x,RDF100max{tind,rind,mind},'r');
%             
%             
%             saveas(gcf,['RDFstdminmax100T' my_num2str(t)...
%                 'rho' my_num2str(r) 'm' num2str(m) '.fig']);
%             saveas(gcf,['RDFstdminmax100T' my_num2str(t)...
%                 'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
%             close all;
%             
%             
%             plot(x,RDF1k{tind,rind,mind});            
%             hold on;
%             plot(x,RDF1k{tind,rind,mind}+...
%                 (RDF1k{tind,rind,mind}.*RDF1kstd{tind,rind,mind})/sqrt(1000),'-r');
%             
%             plot(x,RDF1k{tind,rind,mind}-...
%                 (RDF1k{tind,rind,mind}.*RDF1kstd{tind,rind,mind})/sqrt(1000),'-r');
%             
%             
%             plot(x,RDF1kmin{tind,rind,mind},'r');
%             plot(x,RDF1kmax{tind,rind,mind},'r');
%             
%             
%             saveas(gcf,['RDFstdminmax1kT' my_num2str(t)...
%                 'rho' my_num2str(r) 'm' num2str(m) '.fig']);
%             saveas(gcf,['RDFstdminmax1kT' my_num2str(t)...
%                 'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
%             close all;
%             
%             
%             plot(x,RDF2k{tind,rind,mind});            
%             hold on;
%             plot(x,RDF2k{tind,rind,mind}+...
%                 (RDF2k{tind,rind,mind}.*RDF2kstd{tind,rind,mind})/sqrt(2000),'-r');
%             
%             plot(x,RDF2k{tind,rind,mind}-...
%                 (RDF2k{tind,rind,mind}.*RDF2kstd{tind,rind,mind})/sqrt(2000),'-r');
%             
%             
%             plot(x,RDF2kmin{tind,rind,mind},'r');
%             plot(x,RDF2kmax{tind,rind,mind},'r');
%             
%             
%             saveas(gcf,['RDFstdminmax2kT' my_num2str(t)...
%                 'rho' my_num2str(r) 'm' num2str(m) '.fig']);
%             saveas(gcf,['RDFstdminmax2kT' my_num2str(t)...
%                 'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
%             close all;
%             
%             
%             plot(x,RDF6k{tind,rind,mind});            
%             hold on;
%             plot(x,RDF6k{tind,rind,mind}+...
%                 (RDF6k{tind,rind,mind}.*RDF6kstd{tind,rind,mind})/sqrt(6000),'-r');
%             
%             plot(x,RDF6k{tind,rind,mind}-...
%                 (RDF6k{tind,rind,mind}.*RDF6kstd{tind,rind,mind})/sqrt(6000),'-r');
%             
%             
%             plot(x,RDF6kmin{tind,rind,mind},'r');
%             plot(x,RDF6kmax{tind,rind,mind},'r');
%             
%             
%             saveas(gcf,['RDFstdminmax6kT' my_num2str(t)...
%                 'rho' my_num2str(r) 'm' num2str(m) '.fig']);
%             saveas(gcf,['RDFstdminmax6kT' my_num2str(t)...
%                 'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
%             close all;
%             
%             
%             close all;
        toc
        end 
    end
end

save('all_rdf_minmax.mat');
