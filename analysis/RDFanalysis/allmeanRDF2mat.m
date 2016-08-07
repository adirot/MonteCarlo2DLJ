steps = {'10';'50'; '100'; '1k'; '2k'; '6k'; '12k'};
%start = 4000;

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
            
            if ~and(and(tind == 1, rind == 1),mind == 1)
            list = dir(['*T' my_num2str(t) 'rho' my_num2str(r) '*_m' num2str(m)...
                '*mat']);
            list = {list.name}';

            if length(list) > 1
                disp('more than one file!');
            end

            M = MC2DLJoutput(list{1,1});

            if and(and(t == 0.45,r == 0.005), m==4)
                x = M.data.RDFbins;
            end
            
            RDF10{tind,rind,mind} = mean(M.data.RDFhisto(1:10,1:300));
            RDF50{tind,rind,mind} = mean(M.data.RDFhisto(1:50,1:300));
            RDF100{tind,rind,mind} = mean(M.data.RDFhisto(1:100,1:300));
            RDF1k{tind,rind,mind} = mean(M.data.RDFhisto(1:1000,1:300));
            RDF2k{tind,rind,mind} = mean(M.data.RDFhisto(1:2000,1:300));
            RDF6k{tind,rind,mind} = mean(M.data.RDFhisto(1:6000,1:300));
            
            logRDF10{tind,rind,mind} = -log(RDF10{tind,rind,mind});
            logRDF50{tind,rind,mind} = -log(RDF50{tind,rind,mind});
            logRDF100{tind,rind,mind} = -log(RDF100{tind,rind,mind});
            logRDF1k{tind,rind,mind} = -log(RDF1k{tind,rind,mind});
            logRDF2k{tind,rind,mind} = -log(RDF2k{tind,rind,mind});
            logRDF6k{tind,rind,mind} = -log(RDF6k{tind,rind,mind});
            try
                RDF12k{tind,rind,mind} = mean(M.data.RDFhisto(1:12000,1:300));
                logRDF12k{tind,rind,mind} = -log(mean(M.data.RDFhisto(1:12000,1:300)));
            catch
                RDF12k{tind,rind,mind} = [];
                logRDF12k{tind,rind,mind} = [];
            end
            
            goodInd = ~or(isnan(logRDF10{tind,rind,mind}), logRDF10{tind,rind,mind} == inf);
            logRDF10{tind,rind,mind} = logRDF10{tind,rind,mind}(goodInd);
            xlogRDF10{tind,rind,mind} = x(goodInd);
        
            
            goodInd = ~or(isnan(logRDF50{tind,rind,mind}), logRDF50{tind,rind,mind} == inf);
            logRDF50{tind,rind,mind} = logRDF50{tind,rind,mind}(goodInd);
            xlogRDF50{tind,rind,mind} = x(goodInd);
            
            
            goodInd = ~or(isnan(logRDF100{tind,rind,mind}), logRDF100{tind,rind,mind} == inf);
            logRDF100{tind,rind,mind} = logRDF100{tind,rind,mind}(goodInd);
            xlogRDF100{tind,rind,mind} = x(goodInd);
            
            
            goodInd = ~or(isnan(logRDF1k{tind,rind,mind}), logRDF1k{tind,rind,mind} == inf);
            logRDF1k{tind,rind,mind} = logRDF1k{tind,rind,mind}(goodInd);
            xlogRDF1k{tind,rind,mind} = x(goodInd);

            
            goodInd = ~or(isnan(logRDF2k{tind,rind,mind}), logRDF2k{tind,rind,mind} == inf);
            logRDF2k{tind,rind,mind} = logRDF2k{tind,rind,mind}(goodInd);
            xlogRDF2k{tind,rind,mind} = x(goodInd);
            
            
            goodInd = ~or(isnan(logRDF6k{tind,rind,mind}), logRDF6k{tind,rind,mind} == inf);
            logRDF6k{tind,rind,mind} = logRDF6k{tind,rind,mind}(goodInd);
            xlogRDF6k{tind,rind,mind} = x(goodInd);

            try
                goodInd = ~or(isnan(logRDF12k{tind,rind,mind}), logRDF12k{tind,rind,mind} == inf);
                logRDF12k{tind,rind,mind} = logRDF12k{tind,rind,mind}(goodInd);
                xlogRDF12k{tind,rind,mind} = x(goodInd);
            catch
                xlogRDF12k{tind,rind,mind} = [];
            end

            end
        end 
    end
end

save('logRDFdata2.mat','logRDF10','xlogRDF10','logRDF50','xlogRDF50','logRDF100','xlogRDF100',...
    'logRDF1k','xlogRDF1k','logRDF2k','xlogRDF2k','logRDF6k','xlogRDF6k',...
    'logRDF12k','xlogRDF12k','RDF10','RDF50','RDF100',...
    'RDF1k','RDF2k','RDF6k',...
    'RDF12k');

