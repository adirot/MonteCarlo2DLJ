% fit all matfiles to the potantail (logRDF) for different steps

steps = {'4010';'4050';'4100';'5k';'10k';'14k';'21k';'24k';'28k';'34k'};
% {'10';'50'; '100'; '1k'; '2k'; '6k'; '12k'};
folderName = 'hcr/';




list = dir([folderName 'N*mat']);
list = {list.name}';
start = 4000;
load('logRDFdataHCR.mat');

tind = 0;
rind = 0;
mind = 0;
first = true;

for t = [0.01,0.02,0.03,0.04,0.05,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,1,1.5,2,5]
%for t = [0.01,0.02,0.03,0.04,0.05,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,1,1.5,5]

    tind = tind + 1;
    rind = 0;
    for r = [0.005,0.01,0.05,0.1,0.2]
    %for r = 0.005
    
        rind = rind + 1;
        mind = 0;
        for m = [3 4 5 6]
            mind = mind + 1; 
            
            disp(t); 
            disp(r);
            disp(m);
            disp('-----');
            
            list = dir([folderName '*T' my_num2str(t)...
                'rho' my_num2str(r) '*_m' num2str(m) '*mat']);
            list = {list.name}';

            if length(list) > 1
                disp('more than one file!');
            end
            
            try
                M = MC2DLJoutput([folderName list{1,1}]);
            
                if first
                    x = M.data.RDFbins;
                    first = false;
                end
                missing = false;
            catch
                disp('missing file!');
                missing = true;
            end
            
            if ~missing
                tic;     
        
%             RDF10 = mean(M.data.RDFhisto(start:(start+10),1:300));
%             RDF50 = mean(M.data.RDFhisto(start:(start+50),1:300));
%             RDF100 = mean(M.data.RDFhisto(start:(start+100),1:300));
%             RDF1k = mean(M.data.RDFhisto(start:(start+1000),1:300));
%             RDF2k = mean(M.data.RDFhisto(start:(start+2000),1:300));
%             RDF6k = mean(M.data.RDFhisto(start:(start+6000),1:300));
%             RDF12k = mean(M.data.RDFhisto(start:(start+12000),1:300));
                %steps = {'4010';'4050';'4100';'5k';'10k';'14k';'21k';'24k';'28k';'34k'};

               ys = [logRDF4010{tind,rind,mind};logRDF4050{tind,rind,mind};...
                   logRDF4100{tind,rind,mind};logRDF5k{tind,rind,mind};...
                   logRDF10k{tind,rind,mind};logRDF14k{tind,rind,mind};...
                   logRDF21k{tind,rind,mind};logRDF24k{tind,rind,mind};...
                    logRDF28k{tind,rind,mind};logRDF34k{tind,rind,mind}];
               xs = [xlogRDF4010{tind,rind,mind};xlogRDF4050{tind,rind,mind};...
                   xlogRDF4100{tind,rind,mind};xlogRDF5k{tind,rind,mind};...
                   xlogRDF10k{tind,rind,mind};xlogRDF14k{tind,rind,mind};...
                   xlogRDF21k{tind,rind,mind};xlogRDF24k{tind,rind,mind};...
                    xlogRDF28k{tind,rind,mind};xlogRDF34k{tind,rind,mind}];
               
                % with n set to 12
                %[fitresult,mfit,merror, gof] = createFitRDF(xs, ys, t, r,m,steps);

                % with free n
                %[fitresult,mfit,merror,nfit,nerror,Tfit,Terror, gof] =...
                %    createFitRDF(xs, ys, t, r,m,steps,'freeTandn',true);

                % hcr
                [fitresult, mfit, merror, nfit, nError, Tfit, Terror, gof] =...
                    createFitRDF(xs, ys, t, r, m, steps, 'hcr',true,...
                    'ignoreZerosAtend',true);
                
                try
                    allfittedm4010{tind,rind,mind} = mfit{1};
                    allgof4010{tind,rind,mind} = gof{1};
                catch 
                    allfittedm4010{tind,rind,mind} = [];
                    allgof4010{tind,rind,mind} = [];
                end
                try
                    allfittedm4050{tind,rind,mind} = mfit{2};
                    allgof4050{tind,rind,mind} = gof{2};
                catch 
                    allfittedm4050{tind,rind,mind} = [];
                    allgof4050{tind,rind,mind} = [];
                end
                try
                    allfittedm4100{tind,rind,mind} = mfit{3};
                    allgof4100{tind,rind,mind} = gof{3};
                catch 
                    allfittedm4100{tind,rind,mind} = [];
                    allgof4100{tind,rind,mind} = [];
                end
                try
                    allfittedm5k{tind,rind,mind} = mfit{4};
                    allgof5k{tind,rind,mind} = gof{4};
                catch 
                    allfittedm5k{tind,rind,mind} = [];
                    allgof5k{tind,rind,mind} = [];
                end
                try
                    allfittedm10k{tind,rind,mind} = mfit{5};
                    allgof10k{tind,rind,mind} = gof{5};
                catch 
                    allfittedm10k{tind,rind,mind} = [];
                    allgof10k{tind,rind,mind} = [];
                end
                
                try 
                    allfittedm14k{tind,rind,mind} = mfit{6};
                    allgof14k{tind,rind,mind} = gof{6};
                catch 
                    allfittedm14k{tind,rind,mind} = [];
                    allgof14k{tind,rind,mind} = [];
                end
                
                try
                    allfittedm21k{tind,rind,mind} = mfit{7};
                    allgof21k{tind,rind,mind} = gof{7};
                catch 
                    allfittedm21k{tind,rind,mind} = [];
                    allgof21k{tind,rind,mind} = [];
                end
                
                try
                    allfittedm24k{tind,rind,mind} = mfit{8};
                    allgof24k{tind,rind,mind} = gof{8};
                catch 
                    allfittedm24k{tind,rind,mind} = [];
                    allgof24k{tind,rind,mind} = [];
                end
                
                try
                    allfittedm28k{tind,rind,mind} = mfit{9};
                    allgof28k{tind,rind,mind} = gof{9};
                catch 
                    allfittedm28k{tind,rind,mind} = [];
                    allgof28k{tind,rind,mind} = [];
                end
                
                try
                    allfittedm34k{tind,rind,mind} = mfit{10};
                    allgof34k{tind,rind,mind} = gof{10};
                catch 
                    allfittedm34k{tind,rind,mind} = [];
                    allgof34k{tind,rind,mind} = [];
                end
                

                try
                    allmerror4010{tind,rind,mind} = merror{1};
                catch 
                    allmerror4010{tind,rind,mind} = [];
                end
                try
                    allmerror4050{tind,rind,mind} = merror{2};
                catch 
                    allmerror4050{tind,rind,mind} = [];
                end
                try
                    allmerror4100{tind,rind,mind} = merror{3};
                catch 
                    allmerror4100{tind,rind,mind} = [];
                end
                try
                    allmerror5k{tind,rind,mind} = merror{4};
                catch 
                    allmerror5k{tind,rind,mind} = [];
                end
                
                try
                    allmerror10k{tind,rind,mind} = merror{5};
                catch 
                    allmerror10k{tind,rind,mind} = [];
                end
                
                try
                    allmerror14k{tind,rind,mind} = merror{6};
                catch 
                    allmerror14k{tind,rind,mind} = [];
                end
                
                try
                    allmerror21k{tind,rind,mind} = merror{7};
                catch 
                    allmerror21k{tind,rind,mind} = [];
                end
                
                try
                    allmerror24k{tind,rind,mind} = merror{8};
                catch 
                    allmerror24k{tind,rind,mind} = [];
                end
                
                try
                    allmerror28k{tind,rind,mind} = merror{9};
                catch 
                    allmerror28k{tind,rind,mind} = [];
                end
                
                try
                    allmerror34k{tind,rind,mind} = merror{10};
                catch 
                    allmerror34k{tind,rind,mind} = [];
                end
                

    %             allfittedn10k{tind,rind,mind} = nfit(1);
    %             allfittedn14k{tind,rind,mind} = nfit(2);
    %             allfittedn21k{tind,rind,mind} = nfit(3);
    %             allfittedn24k{tind,rind,mind} = nfit(4);
    %             allfittedn28k{tind,rind,mind} = nfit(5);
    %             allfittedn34k{tind,rind,mind} = nfit(6);
    %             
    %             allnerror10k{tind,rind,mind} = nerror(1);
    %             allnerror14k{tind,rind,mind} = nerror(2);
    %             allnerror21k{tind,rind,mind} = nerror(3);
    %             allnerror24k{tind,rind,mind} = nerror(4);
    %             allnerror28k{tind,rind,mind} = nerror(5);
    %             allnerror34k{tind,rind,mind} = nerror(6);
    %             allnerror128k{tind,rind,mind} = nerror(7);

                
                try
                    allfittedT4010{tind,rind,mind} = Tfit{1};
                catch 
                    allfittedT4010{tind,rind,mind} = [];
                end
                 try
                    allfittedT4050{tind,rind,mind} = Tfit{2};
                catch 
                    allfittedT4050{tind,rind,mind} = [];
                 end
                try
                    allfittedT4100{tind,rind,mind} = Tfit{3};
                catch 
                    allfittedT4100{tind,rind,mind} = [];
                end
                 
                 try
                    allfittedT5k{tind,rind,mind} = Tfit{4};
                catch 
                    allfittedT5k{tind,rind,mind} = [];
                end
                 
                try
                    allfittedT10k{tind,rind,mind} = Tfit{5};
                catch 
                    allfittedT10k{tind,rind,mind} = [];
                end
                    
                try
                    allfittedT14k{tind,rind,mind} = Tfit{6};
                catch 
                    allfittedT14k{tind,rind,mind} = [];
                end
                    
                try
                    allfittedT21k{tind,rind,mind} = Tfit{7};
                catch 
                    allfittedT21k{tind,rind,mind} = [];
                end
                
                
                try
                    allfittedT24k{tind,rind,mind} = Tfit{8};
                catch 
                    allfittedT24k{tind,rind,mind} = [];
                end
                
                
                try    
                    allfittedT28k{tind,rind,mind} = Tfit{9};
                catch 
                    allfittedT28k{tind,rind,mind} = [];
                end
                
                
                try   
                    allfittedT34k{tind,rind,mind} = Tfit{10};
                catch 
                    allfittedT34k{tind,rind,mind} = [];
                end
                try    
                    allTerror4010{tind,rind,mind} = Terror(1,1:2);
                catch 
                    allTerror4010{tind,rind,mind} = [];
                end
                try    
                    allTerror4050{tind,rind,mind} = Terror(2,1:2);
                catch 
                    allTerror4050{tind,rind,mind} = [];
                end
                try    
                    allTerror4100{tind,rind,mind} = Terror(3,1:2);
                catch 
                    allTerror4100{tind,rind,mind} = [];
                end
                try    
                    allTerror5k{tind,rind,mind} = Terror(4,1:2);
                catch 
                    allTerror5k{tind,rind,mind} = [];
                end
                
                try    
                    allTerror10k{tind,rind,mind} = Terror(5,1:2);
                catch 
                    allTerror10k{tind,rind,mind} = [];
                end
                
                try
                    allTerror14k{tind,rind,mind} = Terror(6,1:2);
                catch 
                    allTerror14k{tind,rind,mind} = [];
                end
                
                try   
                    allTerror21k{tind,rind,mind} = Terror(7,1:2);
                catch 
                    allTerror21k{tind,rind,mind} = [];
                end
                
                try   
                    allTerror24k{tind,rind,mind} = Terror(8,1:2);
                catch 
                    allTerror24k{tind,rind,mind} = [];
                end
                
                try   
                    allTerror28k{tind,rind,mind} = Terror(9,1:2);
                catch 
                    allTerror28k{tind,rind,mind} = [];
                end
                
                try   
                    allTerror34k{tind,rind,mind} = Terror(10,1:2);
                catch 
                    allTerror34k{tind,rind,mind} = [];
                end
                

                saveas(gcf,[folderName 'fitRDFT' my_num2str(t)...
                    'rho' my_num2str(r) 'm' num2str(m) 'freenT_steps.fig']);
                saveas(gcf,[folderName 'fitRDFT' my_num2str(t)...
                    'rho' my_num2str(r) 'm' num2str(m) 'freenT_steps.jpg']);
                close all;

                toc;
            end
        end 
    end
end

save([folderName 'all_fittedm_Thcr.mat']);
