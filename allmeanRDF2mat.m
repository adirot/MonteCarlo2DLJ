

steps = {'10';'50'; '100'; '1k'; '2k'; '6k'; '12k'};
list = dir('N*mat');
list = {list.name}';
start = 4000;

tind = 0;
rind = 0;
mind = 0;
nind = 0;

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
            
        
            logRDF10{tind,rind,mind} = -log(mean(M.data.RDFhisto(start:(start+10),1:300)));
            logRDF50{tind,rind,mind} = -log(mean(M.data.RDFhisto(start:(start+50),1:300)));
            logRDF100{tind,rind,mind} = -log(mean(M.data.RDFhisto(start:(start+100),1:300)));
            logRDF1k{tind,rind,mind} = -log(mean(M.data.RDFhisto(start:(start+1000),1:300)));
            logRDF2k{tind,rind,mind} = -log(mean(M.data.RDFhisto(start:(start+2000),1:300)));
            logRDF6k{tind,rind,mind} = -log(mean(M.data.RDFhisto(start:(start+6000),1:300)));
            logRDF12k{tind,rind,mind} = -log(mean(M.data.RDFhisto(start:(start+12000),1:300)));
        

        end 
    end
end

save('logRDFdata.mat');

