steps = {'4010';'4050';'4100';'5k';'10k';'14k';'21k';'24k';'28k';'34k'};
start = 4000;
folderName = 'hcr/'; 
x = (1:300)*10/300;

tind = 0;
rind = 0;
mind = 0;

for t = [0.01,0.02,0.03,0.04,0.05,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,1,1.5,2,5]
%for t = [0.05,0.07,0.1,0.2,0.3,0.4,2]
%for t = [0.08 0.09 0.3]
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
            
            list = dir([folderName '*T' my_num2str(t) 'rho'...
                my_num2str(r) '*_m' num2str(m) '*mat']);
            list = {list.name}';

            if length(list) > 1
                disp('more than one file!');
                
            end
            
            if isempty(list)
                disp('missing file!');
            else
                M = [];
                listind = 1; maxindIndata = 10;
                while  listind <= length(list)
                    try 
                        M1 = MC2DLJoutput([folderName list{listind,1}]);
                        if M1.indIndata > maxindIndata
                            M = M1;
                            maxindIndata = M1.indIndata;
                        end
                        listind = listind + 1;
                    catch
                        listind = listind + 1;
                    end
                end
                
                if ~isempty(M)
                    
                    try
                        RDF4010{tind,rind,mind} = mean(M.data.RDFhisto(1,1:300,start:4010),3);
                        logRDF4010{tind,rind,mind} = -log(RDF4010{tind,rind,mind});
                    catch
                        RDF4010{tind,rind,mind} = zeros(1,300);
                        logRDF4010{tind,rind,mind} = zeros(1,300);
                    end
                    try
                        RDF4050{tind,rind,mind} = mean(M.data.RDFhisto(1,1:300,start:4050),3);
                        logRDF4050{tind,rind,mind} = -log(RDF4050{tind,rind,mind});
                    catch
                        RDF4050{tind,rind,mind} = zeros(1,300);
                        logRDF4050{tind,rind,mind} = zeros(1,300);
                    end


                    try
                        RDF4100{tind,rind,mind} = mean(M.data.RDFhisto(1,1:300,start:4100),3);
                        logRDF4100{tind,rind,mind} = -log(RDF4100{tind,rind,mind});
                    catch
                        RDF4100{tind,rind,mind} = zeros(1,300);
                        logRDF4100{tind,rind,mind} = zeros(1,300);
                    end


                    try
                        RDF5k{tind,rind,mind} = mean(M.data.RDFhisto(1,1:300,start:5000),3);
                        logRDF5k{tind,rind,mind} = -log(RDF5k{tind,rind,mind});
                    catch
                        RDF5k{tind,rind,mind} = zeros(1,300);
                        logRDF5k{tind,rind,mind} = zeros(1,300);
                    end

                    try
                        RDF10k{tind,rind,mind} = mean(M.data.RDFhisto(1,1:300,start:10000),3);
                        logRDF10k{tind,rind,mind} = -log(RDF10k{tind,rind,mind});
                    catch
                        RDF10k{tind,rind,mind} = zeros(1,300);
                        logRDF10k{tind,rind,mind} = zeros(1,300);
                    end

                    try
                        RDF14k{tind,rind,mind} = mean(M.data.RDFhisto(1,1:300,start:14000),3);
                        logRDF14k{tind,rind,mind} = -log(RDF14k{tind,rind,mind});
                    catch
                        RDF14k{tind,rind,mind} = zeros(1,300);
                        logRDF14k{tind,rind,mind} = zeros(1,300);

                    end

                    try
                        RDF21k{tind,rind,mind} = mean(M.data.RDFhisto(1,1:300,start:21000),3);
                        logRDF21k{tind,rind,mind} = -log(RDF21k{tind,rind,mind});

                    catch
                        RDF21k{tind,rind,mind} = zeros(1,300);
                        logRDF21k{tind,rind,mind} = zeros(1,300);

                    end

                    try
                        RDF24k{tind,rind,mind} = mean(M.data.RDFhisto(1,1:300,start:24000),3);
                        logRDF24k{tind,rind,mind} = -log(RDF24k{tind,rind,mind});

                    catch
                        RDF24k{tind,rind,mind} = zeros(1,300);
                        logRDF24k{tind,rind,mind} = zeros(1,300);

                    end

                    try
                        RDF28k{tind,rind,mind} = mean(M.data.RDFhisto(1,1:300,start:28000),3);
                        logRDF28k{tind,rind,mind} = -log(RDF28k{tind,rind,mind});

                    catch
                        RDF28k{tind,rind,mind} = zeros(1,300);
                        logRDF28k{tind,rind,mind} = zeros(1,300);

                    end

                    try
                        RDF34k{tind,rind,mind} = mean(M.data.RDFhisto(1,1:300,start:34000),3);
                        logRDF34k{tind,rind,mind} = -log(RDF34k{tind,rind,mind});

                    catch
                        RDF34k{tind,rind,mind} = zeros(1,300);
                        logRDF34k{tind,rind,mind} = zeros(1,300);
                    end

                    try

                        goodInd = ~or(isnan(logRDF4010{tind,rind,mind}), logRDF4010{tind,rind,mind} == inf);
                        logRDF4010{tind,rind,mind} = logRDF4010{tind,rind,mind}(goodInd);
                        xlogRDF4010{tind,rind,mind} = x(goodInd);
                        xlogRDF4010{tind,rind,mind} = [xlogRDF4010{tind,rind,mind} zeros(1,300-length(xlogRDF4010{tind,rind,mind}))];
                        logRDF4010{tind,rind,mind} = [logRDF4010{tind,rind,mind} zeros(1,300-length(logRDF4010{tind,rind,mind}))] ;
                    catch
                        logRDF4010{tind,rind,mind} = zeros(1,300);
                        xlogRDF4010{tind,rind,mind} = zeros(1,300);
                    end
                    try

                        goodInd = ~or(isnan(logRDF4050{tind,rind,mind}), logRDF4050{tind,rind,mind} == inf);
                        logRDF4050{tind,rind,mind} = logRDF4050{tind,rind,mind}(goodInd);
                        xlogRDF4050{tind,rind,mind} = x(goodInd);
                        xlogRDF4050{tind,rind,mind} = [xlogRDF4050{tind,rind,mind} zeros(1,300-length(xlogRDF4050{tind,rind,mind}))];
                        logRDF4050{tind,rind,mind} = [logRDF4050{tind,rind,mind} zeros(1,300-length(logRDF4050{tind,rind,mind}))] ;
                    catch
                        logRDF4050{tind,rind,mind} = zeros(1,300);
                        xlogRDF4050{tind,rind,mind} = zeros(1,300);
                    end
                    try

                        goodInd = ~or(isnan(logRDF4100{tind,rind,mind}), logRDF4100{tind,rind,mind} == inf);
                        logRDF4100{tind,rind,mind} = logRDF4100{tind,rind,mind}(goodInd);
                        xlogRDF4100{tind,rind,mind} = x(goodInd);
                        xlogRDF4100{tind,rind,mind} = [xlogRDF4100{tind,rind,mind} zeros(1,300-length(xlogRDF4100{tind,rind,mind}))];
                        logRDF4100{tind,rind,mind} = [logRDF4100{tind,rind,mind} zeros(1,300-length(logRDF4100{tind,rind,mind}))] ;
                    catch
                        logRDF4100{tind,rind,mind} = zeros(1,300);
                        xlogRDF4100{tind,rind,mind} = zeros(1,300);
                    end
                    try

                        goodInd = ~or(isnan(logRDF5k{tind,rind,mind}), logRDF5k{tind,rind,mind} == inf);
                        logRDF5k{tind,rind,mind} = logRDF5k{tind,rind,mind}(goodInd);
                        xlogRDF5k{tind,rind,mind} = x(goodInd);
                        xlogRDF5k{tind,rind,mind} = [xlogRDF5k{tind,rind,mind} zeros(1,300-length(xlogRDF5k{tind,rind,mind}))];
                        logRDF5k{tind,rind,mind} = [logRDF5k{tind,rind,mind} zeros(1,300-length(logRDF5k{tind,rind,mind}))] ;
                    catch
                        logRDF5k{tind,rind,mind} = zeros(1,300);
                        xlogRDF5k{tind,rind,mind} = zeros(1,300);
                    end

                    try

                        goodInd = ~or(isnan(logRDF10k{tind,rind,mind}), logRDF10k{tind,rind,mind} == inf);
                        logRDF10k{tind,rind,mind} = logRDF10k{tind,rind,mind}(goodInd);
                        xlogRDF10k{tind,rind,mind} = x(goodInd);
                        xlogRDF10k{tind,rind,mind} = [xlogRDF10k{tind,rind,mind} zeros(1,300-length(xlogRDF10k{tind,rind,mind}))];
                        logRDF10k{tind,rind,mind} = [logRDF10k{tind,rind,mind} zeros(1,300-length(logRDF10k{tind,rind,mind}))] ;
                    catch
                        logRDF10k{tind,rind,mind} = zeros(1,300);
                        xlogRDF10k{tind,rind,mind} = zeros(1,300);
                    end
                    try
                        goodInd = ~or(isnan(logRDF14k{tind,rind,mind}), logRDF14k{tind,rind,mind} == inf);
                        logRDF14k{tind,rind,mind} = logRDF14k{tind,rind,mind}(goodInd);
                        xlogRDF14k{tind,rind,mind} = x(goodInd);
                        xlogRDF14k{tind,rind,mind} = [xlogRDF14k{tind,rind,mind} zeros(1,300-length(xlogRDF14k{tind,rind,mind}))];
                        logRDF14k{tind,rind,mind} = [logRDF14k{tind,rind,mind} zeros(1,300-length(logRDF14k{tind,rind,mind}))] ;
                    catch
                        logRDF14k{tind,rind,mind} = zeros(1,300);
                        xlogRDF14k{tind,rind,mind} = zeros(1,300);
                    end

                    try
                        goodInd = ~or(isnan(logRDF21k{tind,rind,mind}), logRDF21k{tind,rind,mind} == inf);
                        logRDF21k{tind,rind,mind} = logRDF21k{tind,rind,mind}(goodInd);
                        xlogRDF21k{tind,rind,mind} = x(goodInd);
                        xlogRDF21k{tind,rind,mind} = [xlogRDF21k{tind,rind,mind} zeros(1,300-length(xlogRDF21k{tind,rind,mind}))];
                        logRDF21k{tind,rind,mind} = [logRDF21k{tind,rind,mind} zeros(1,300-length(logRDF21k{tind,rind,mind}))] ;
                    catch
                        logRDF21k{tind,rind,mind} = zeros(1,300);
                        xlogRDF21k{tind,rind,mind} = zeros(1,300);
                    end


                    try
                        goodInd = ~or(isnan(logRDF24k{tind,rind,mind}), logRDF24k{tind,rind,mind} == inf);
                        logRDF24k{tind,rind,mind} = logRDF24k{tind,rind,mind}(goodInd);
                        xlogRDF24k{tind,rind,mind} = x(goodInd);
                        xlogRDF24k{tind,rind,mind} = [xlogRDF24k{tind,rind,mind} zeros(1,300-length(xlogRDF24k{tind,rind,mind}))];
                        logRDF24k{tind,rind,mind} = [logRDF24k{tind,rind,mind} zeros(1,300-length(logRDF24k{tind,rind,mind}))] ;
                    catch
                        logRDF24k{tind,rind,mind} = zeros(1,300);
                        xlogRDF24k{tind,rind,mind} = zeros(1,300);
                    end

                    try
                        goodInd = ~or(isnan(logRDF28k{tind,rind,mind}), logRDF28k{tind,rind,mind} == inf);
                        logRDF28k{tind,rind,mind} = logRDF28k{tind,rind,mind}(goodInd);
                        xlogRDF28k{tind,rind,mind} = x(goodInd);
                        xlogRDF28k{tind,rind,mind} = [xlogRDF28k{tind,rind,mind} zeros(1,300-length(xlogRDF28k{tind,rind,mind}))];
                        logRDF28k{tind,rind,mind} = [logRDF28k{tind,rind,mind} zeros(1,300-length(logRDF28k{tind,rind,mind}))] ;
                    catch
                        logRDF28k{tind,rind,mind} = zeros(1,300);
                        xlogRDF28k{tind,rind,mind} = zeros(1,300);
                    end
                    try
                        goodInd = ~or(isnan(logRDF34k{tind,rind,mind}), logRDF34k{tind,rind,mind} == inf);
                        logRDF34k{tind,rind,mind} = logRDF34k{tind,rind,mind}(goodInd);
                        xlogRDF34k{tind,rind,mind} = x(goodInd);
                        xlogRDF34k{tind,rind,mind} = [xlogRDF34k{tind,rind,mind} zeros(1,300-length(xlogRDF34k{tind,rind,mind}))];
                        logRDF34k{tind,rind,mind} = [logRDF34k{tind,rind,mind} zeros(1,300-length(logRDF34k{tind,rind,mind}))] ;
                    catch
                        logRDF34k{tind,rind,mind} = zeros(1,300);
                        xlogRDF34k{tind,rind,mind} = zeros(1,300);
                    end
                end

            end
        end 
    end
end

save('logRDFdataHCR.mat','logRDF4010','xlogRDF4010','logRDF4050',...
    'xlogRDF4050','logRDF4100','xlogRDF4100',...
    'logRDF5k','xlogRDF5k','logRDF10k','xlogRDF10k','logRDF14k',...
    'xlogRDF14k','logRDF21k','xlogRDF21k',...
    'logRDF24k','xlogRDF24k','logRDF28k','xlogRDF28k',...
    'logRDF34k','xlogRDF34k',...
    'RDF4010','RDF4050','RDF4100','RDF5k','RDF10k','RDF14k','RDF21k',...
    'RDF24k','RDF28k','RDF34k');

