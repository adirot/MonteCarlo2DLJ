%% get std of -log(RDF) %%

load('all_rdf_minmax.mat');
load('logRDFdata2.mat');
load('all_multifitted.mat');
steps = {'10';'50'; '100'; '1k'; '2k'; '6k'; '12k'};

tind = 0;

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
            
            if ~and(and(tind==1,rind==1),mind==1)
            list = dir(['*T' my_num2str(t) 'rho' my_num2str(r) '*_m' num2str(m)...
                '*mat']);
            list = {list.name}';

            if length(list) > 1
                disp('more than one file!');
            end

            M = MC2DLJoutput(list{1,1});

                x = M.data.RDFbins;
            
                 
        
            RDF10std{tind,rind,mind} = std(logRDF10{tind,rind,mind});
            RDF50std{tind,rind,mind} = std(logRDF50{tind,rind,mind});
            RDF100std{tind,rind,mind} = std(logRDF100{tind,rind,mind});
            RDF1kstd{tind,rind,mind} = std(logRDF1k{tind,rind,mind});
            RDF2kstd{tind,rind,mind} = std(logRDF2k{tind,rind,mind});
            RDF6kstd{tind,rind,mind} = std(logRDF6k{tind,rind,mind});
            RDF12kstd{tind,rind,mind} = std(logRDF12k{tind,rind,mind});
        
            
            plot(x,RDF12k{tind,rind,mind});            
            hold on;
            plot(x,RDF12k{tind,rind,mind}+...
                (RDF12k{tind,rind,mind}.*RDF12kstd{tind,rind,mind})/sqrt(12000),'-r');
            
            plot(x,RDF12k{tind,rind,mind}-...
                (RDF12k{tind,rind,mind}.*RDF12kstd{tind,rind,mind})/sqrt(12000),'-r');
            
            plot(x,RDF12kstd{tind,rind,mind},'r');
            plot(x,RDF12kmax{tind,rind,mind},'r');
            
            saveas(gcf,['RDFstdminmax12kT' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.fig']);
            saveas(gcf,['RDFstdminmax12kT' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
            close all;

            
            plot(x,RDF10{tind,rind,mind});            
            hold on;
            plot(x,RDF10{tind,rind,mind}+...
                (RDF10{tind,rind,mind}.*RDF10std{tind,rind,mind})/sqrt(10),'-r');
            
            plot(x,RDF10{tind,rind,mind}-...
                (RDF10{tind,rind,mind}.*RDF10std{tind,rind,mind})/sqrt(10),'-r');
            
            plot(x,RDF10std{tind,rind,mind},'r');
            plot(x,RDF10max{tind,rind,mind},'r');
            
            saveas(gcf,['RDFstdminmax10T' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.fig']);
            saveas(gcf,['RDFstdminmax10T' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
            close all;
            
            
            plot(x,RDF50{tind,rind,mind});            
            hold on;
            plot(x,RDF50{tind,rind,mind}+...
                (RDF50{tind,rind,mind}.*RDF50std{tind,rind,mind})/sqrt(50),'-r');
            
            plot(x,RDF50{tind,rind,mind}-...
                (RDF50{tind,rind,mind}.*RDF50std{tind,rind,mind})/sqrt(50),'-r');
            
            
            plot(x,RDF50std{tind,rind,mind},'r');
            plot(x,RDF50max{tind,rind,mind},'r');
            
            saveas(gcf,['RDFstdminmax50T' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.fig']);
            saveas(gcf,['RDFstdminmax50T' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
            close all;
            
            
            plot(x,RDF100{tind,rind,mind});            
            hold on;
            plot(x,RDF100{tind,rind,mind}+...
                (RDF100{tind,rind,mind}.*RDF100std{tind,rind,mind})/sqrt(100),'-r');
            
            plot(x,RDF100{tind,rind,mind}-...
                (RDF100{tind,rind,mind}.*RDF100std{tind,rind,mind})/sqrt(100),'-r');
            
            
            plot(x,RDF100std{tind,rind,mind},'r');
            plot(x,RDF100max{tind,rind,mind},'r');
            
            
            saveas(gcf,['RDFstdminmax100T' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.fig']);
            saveas(gcf,['RDFstdminmax100T' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
            close all;
            
            
            plot(x,RDF1k{tind,rind,mind});            
            hold on;
            plot(x,RDF1k{tind,rind,mind}+...
                (RDF1k{tind,rind,mind}.*RDF1kstd{tind,rind,mind})/sqrt(1000),'-r');
            
            plot(x,RDF1k{tind,rind,mind}-...
                (RDF1k{tind,rind,mind}.*RDF1kstd{tind,rind,mind})/sqrt(1000),'-r');
            
            
            plot(x,RDF1kstd{tind,rind,mind},'r');
            plot(x,RDF1kmax{tind,rind,mind},'r');
            
            
            saveas(gcf,['RDFstdminmax1kT' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.fig']);
            saveas(gcf,['RDFstdminmax1kT' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
            close all;
            
            
            plot(x,RDF2k{tind,rind,mind});            
            hold on;
            plot(x,RDF2k{tind,rind,mind}+...
                (RDF2k{tind,rind,mind}.*RDF2kstd{tind,rind,mind})/sqrt(2000),'-r');
            
            plot(x,RDF2k{tind,rind,mind}-...
                (RDF2k{tind,rind,mind}.*RDF2kstd{tind,rind,mind})/sqrt(2000),'-r');
            
            
            plot(x,RDF2kstd{tind,rind,mind},'r');
            plot(x,RDF2kmax{tind,rind,mind},'r');
            
            
            saveas(gcf,['RDFstdminmax2kT' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.fig']);
            saveas(gcf,['RDFstdminmax2kT' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
            close all;
            
            
            plot(x,RDF6k{tind,rind,mind});            
            hold on;
            plot(x,RDF6k{tind,rind,mind}+...
                (RDF6k{tind,rind,mind}.*RDF6kstd{tind,rind,mind})/sqrt(6000),'-r');
            
            plot(x,RDF6k{tind,rind,mind}-...
                (RDF6k{tind,rind,mind}.*RDF6kstd{tind,rind,mind})/sqrt(6000),'-r');
            
            
            plot(x,RDF6kstd{tind,rind,mind},'r');
            plot(x,RDF6kmax{tind,rind,mind},'r');
            
            
            saveas(gcf,['RDFstdminmax6kT' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.fig']);
            saveas(gcf,['RDFstdminmax6kT' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
            close all;
            
            
            close all;
            end
        end 
    end
end

save('all_rdf_std.mat');
