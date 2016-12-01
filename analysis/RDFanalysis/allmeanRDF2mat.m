steps = {'10k';'14k';'21k';'24k';'28k';'34k'};
start = 4000;
folderName = 'hcr/'; 

tind = 0;
rind = 0;
mind = 0;

for t = [0.01,0.02,0.03,0.04,0.05,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,1,1.5,2,5]
    tind = tind + 1;
    rind = 0;
    for r = [0.005,0.01,0.05,0.1,0.2]
        rind = rind + 1;
        mind = 0;
        for m = 3
            mind = mind + 1; 
            disp(t); 
            disp(r);
            disp(m);
            disp('-----');
            
            list = dir([folderName '*T' my_num2str(t) 'rho'...
                my_num2str(r) '*_m' num2str(m) '*mat']);
            list = {list.name}';

            if length(list) > 1
                disp('more than one file!');-
                
            end
            
            if isempty(list)
                disp('missing file!');
            else
                M = MC2DLJoutput([folderName list{1,1}]);
                if length(list) > 1
                    if M.indIndata < 1000
                        M = MC2DLJoutput([folderName list{2,1}]);
                    end
                end
            
                if and(and(tind == 1,rind == 1), mind==2)
                    x = M.data.RDFbins;
                end

                try
                    RDF10k{tind,rind,mind} = mean(M.data.RDFhisto(1,1:300,start:10000),3);
                    logRDF10k{tind,rind,mind} = -log(RDF10k{tind,rind,mind});
                catch
                    RDF10k{tind,rind,mind} = [];
                    logRDF10k{tind,rind,mind} = [];
                end

                try
                    RDF14k{tind,rind,mind} = mean(M.data.RDFhisto(1,1:300,start:14000),3);
                    logRDF14k{tind,rind,mind} = -log(RDF14k{tind,rind,mind});
                catch
                    RDF14k{tind,rind,mind} = [];
                    logRDF14k{tind,rind,mind} = [];

                end

                try
                    RDF21k{tind,rind,mind} = mean(M.data.RDFhisto(1,1:300,start:21000),3);
                    logRDF21k{tind,rind,mind} = -log(RDF21k{tind,rind,mind});

                catch
                    RDF21k{tind,rind,mind} = [];
                    logRDF21k{tind,rind,mind} = [];

                end

                try
                    RDF24k{tind,rind,mind} = mean(M.data.RDFhisto(1,1:300,start:24000),3);
                    logRDF24k{tind,rind,mind} = -log(RDF24k{tind,rind,mind});

                catch
                    RDF24k{tind,rind,mind} = [];
                    logRDF24k{tind,rind,mind} = [];

                end

                try
                    RDF28k{tind,rind,mind} = mean(M.data.RDFhisto(1,1:300,start:28000),3);
                    logRDF28k{tind,rind,mind} = -log(RDF28k{tind,rind,mind});

                catch
                    RDF28k{tind,rind,mind} = [];
                    logRDF28k{tind,rind,mind} = [];

                end

                try
                    RDF34k{tind,rind,mind} = mean(M.data.RDFhisto(1,1:300,start:34000),3);
                    logRDF34k{tind,rind,mind} = -log(RDF34k{tind,rind,mind});

                catch
                    RDF34k{tind,rind,mind} = [];
                    logRDF34k{tind,rind,mind} = [];
                end

                %'10k';'14k';'21k';'24k';'28k';'34k'
                try
                    goodInd = ~or(isnan(logRDF10k{tind,rind,mind}), logRDF10k{tind,rind,mind} == inf);
                    logRDF10k{tind,rind,mind} = logRDF10k{tind,rind,mind}(goodInd);
                    xlogRDF10k{tind,rind,mind} = x(goodInd);
                catch
                    xlogRDF10k{tind,rind,mind} = [];
                end

                try
                    goodInd = ~or(isnan(logRDF14k{tind,rind,mind}), logRDF14k{tind,rind,mind} == inf);
                    logRDF14k{tind,rind,mind} = logRDF14k{tind,rind,mind}(goodInd);
                    xlogRDF14k{tind,rind,mind} = x(goodInd);
                catch
                    xlogRDF14k{tind,rind,mind} = [];
                end

                try
                    goodInd = ~or(isnan(logRDF21k{tind,rind,mind}), logRDF21k{tind,rind,mind} == inf);
                    logRDF21k{tind,rind,mind} = logRDF21k{tind,rind,mind}(goodInd);
                    xlogRDF21k{tind,rind,mind} = x(goodInd);
                catch
                    xlogRDF21k{tind,rind,mind} = [];
                end


                try
                    goodInd = ~or(isnan(logRDF24k{tind,rind,mind}), logRDF24k{tind,rind,mind} == inf);
                    logRDF24k{tind,rind,mind} = logRDF24k{tind,rind,mind}(goodInd);
                    xlogRDF24k{tind,rind,mind} = x(goodInd);
                catch
                    xlogRDF24k{tind,rind,mind} = [];
                end

                try
                    goodInd = ~or(isnan(logRDF28k{tind,rind,mind}), logRDF28k{tind,rind,mind} == inf);
                    logRDF28k{tind,rind,mind} = logRDF28k{tind,rind,mind}(goodInd);
                    xlogRDF28k{tind,rind,mind} = x(goodInd);
                catch
                    xlogRDF28k{tind,rind,mind} = [];
                end

                try
                    goodInd = ~or(isnan(logRDF34k{tind,rind,mind}), logRDF34k{tind,rind,mind} == inf);
                    logRDF34k{tind,rind,mind} = logRDF34k{tind,rind,mind}(goodInd);
                    xlogRDF34k{tind,rind,mind} = x(goodInd);
                catch
                    xlogRDF34k{tind,rind,mind} = [];
                end

            end
        end 
    end
end

save('logRDFdataHCR.mat','logRDF10k','xlogRDF10k','logRDF14k',...
    'xlogRDF14k','logRDF21k','xlogRDF21k',...
    'logRDF24k','xlogRDF24k','logRDF28k','xlogRDF28k',...
    'logRDF34k','xlogRDF34k',...
    'RDF10k','RDF14k','RDF21k',...
    'RDF24k','RDF28k','RDF34k');

