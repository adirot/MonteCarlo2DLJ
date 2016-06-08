%% get all RDFs, not average. get the highst and lowest RDF, %%
%% for extraction of the temperature %%

% fit all matfiles to the potantail (logRDF) for different steps

% list = dir('N*mat');
% list = {list.name}';
% start = 4000;
% 

load('all_RDF_with_std.mat');

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
            
%             list = dir(['*T' my_num2str(t) 'rho' my_num2str(r) '*_m' num2str(m)...
%                 '*mat']);
%             list = {list.name}';

%             if length(list) > 1
%                 disp('more than one file!');
%             end

%             M = MC2DLJoutput(list{1,1});

%             if and(and(t == 0.45,r == 0.005), m==3)
%                 x = M.data.RDFbins;
%             end
            
            tic;     
        
%             RDF10{tind,rind,mind} = mean(M.data.RDFhisto(start:(start+10),1:300));
%             RDF50{tind,rind,mind} = mean(M.data.RDFhisto(start:(start+50),1:300));
%             RDF100{tind,rind,mind} = mean(M.data.RDFhisto(start:(start+100),1:300));
%             RDF1k{tind,rind,mind} = mean(M.data.RDFhisto(start:(start+1000),1:300));
%             RDF2k{tind,rind,mind} = mean(M.data.RDFhisto(start:(start+2000),1:300));
%             RDF6k{tind,rind,mind} = mean(M.data.RDFhisto(start:(start+6000),1:300));
%             RDF12k{tind,rind,mind} = mean(M.data.RDFhisto(start:(start+12000),1:300));
% 
%             RDF10std{tind,rind,mind} = std(M.data.RDFhisto(start:(start+10),1:300));
%             RDF50std{tind,rind,mind} = std(M.data.RDFhisto(start:(start+50),1:300));
%             RDF100std{tind,rind,mind} = std(M.data.RDFhisto(start:(start+100),1:300));
%             RDF1kstd{tind,rind,mind} = std(M.data.RDFhisto(start:(start+1000),1:300));
%             RDF2kstd{tind,rind,mind} = std(M.data.RDFhisto(start:(start+2000),1:300));
%             RDF6kstd{tind,rind,mind} = std(M.data.RDFhisto(start:(start+6000),1:300));
%             RDF12kstd{tind,rind,mind} = std(M.data.RDFhisto(start:(start+12000),1:300));

            plot(x,RDF12k{tind,rind,mind});            
            hold on;
            plot(x,RDF12k{tind,rind,mind}+RDF12kstd{tind,rind,mind},'k');
            plot(x,RDF12k{tind,rind,mind}-RDF12kstd{tind,rind,mind},'k');
            
            saveas(gcf,['RDFstdT' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.fig']);
            saveas(gcf,['RDFstdT' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
            close all;

            
            plot(x,RDF10{tind,rind,mind});            
            hold on;
            plot(x,RDF10{tind,rind,mind}+RDF10std{tind,rind,mind},'k');
            plot(x,RDF10{tind,rind,mind}-RDF10std{tind,rind,mind},'k');
            
            saveas(gcf,['RDFstd10T' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.fig']);
            saveas(gcf,['RDFstd10T' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
            close all;
            
            
            plot(x,RDF50{tind,rind,mind});            
            hold on;
            plot(x,RDF50{tind,rind,mind}+RDF50std{tind,rind,mind},'k');
            plot(x,RDF50{tind,rind,mind}-RDF50std{tind,rind,mind},'k');
            
            
            saveas(gcf,['RDFstd50T' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.fig']);
            saveas(gcf,['RDFstd50T' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
            close all;
            
            
            plot(x,RDF100{tind,rind,mind});            
            hold on;
            plot(x,RDF100{tind,rind,mind}+RDF100std{tind,rind,mind},'k');
            plot(x,RDF100{tind,rind,mind}-RDF100std{tind,rind,mind},'k');
            
            
            saveas(gcf,['RDFstd100T' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.fig']);
            saveas(gcf,['RDFstd100T' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
            close all;
            
            
            plot(x,RDF1k{tind,rind,mind});            
            hold on;
            plot(x,RDF1k{tind,rind,mind}+RDF1kstd{tind,rind,mind},'k');
            plot(x,RDF1k{tind,rind,mind}-RDF1kstd{tind,rind,mind},'k');
            
            
            saveas(gcf,['RDFstd1kT' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.fig']);
            saveas(gcf,['RDFstd1kT' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
            close all;
            
            
            plot(x,RDF2k{tind,rind,mind});            
            hold on;
            
            plot(x,RDF2k{tind,rind,mind}+RDF2kstd{tind,rind,mind},'k');
            plot(x,RDF2k{tind,rind,mind}-RDF2kstd{tind,rind,mind},'k');
            
            saveas(gcf,['RDFstd2kT' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.fig']);
            saveas(gcf,['RDFstd2kT' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
            close all;
            
            
            plot(x,RDF6k{tind,rind,mind});            
            hold on;
            
            plot(x,RDF6k{tind,rind,mind}+RDF6kstd{tind,rind,mind},'k');
            plot(x,RDF6k{tind,rind,mind}-RDF6kstd{tind,rind,mind},'k');
            
            saveas(gcf,['RDFstd6kT' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.fig']);
            saveas(gcf,['RDFstd6kT' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
            close all;
            
            
            plot(x,RDF12k{tind,rind,mind});            
            hold on;
            
            plot(x,RDF12k{tind,rind,mind}+RDF12kstd{tind,rind,mind},'k');
            plot(x,RDF12k{tind,rind,mind}-RDF12kstd{tind,rind,mind},'k');
            
            saveas(gcf,['RDFstd12kT' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.fig']);
            saveas(gcf,['RDFstd12kT' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
            close all;
            
            toc
        end 
    end
end
