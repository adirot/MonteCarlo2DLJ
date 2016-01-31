function percent_done = report(varargin)

        p = inputParser();
        addOptional(p, 'N', ''); 
        addOptional(p, 'fileList', []);
        addOptional(p, 'existsInMat', []);
        addOptional(p, 'monte', false);        
        parse(p, varargin{:});
        Results = p.Results;
        N = Results.N;
        fileList = Results.fileList;
        existsInMat = Results.existsInMat;
        monte = Results.monte;
        
        if isempty(fileList)
            if ~strcmp(N,'')
                N = num2str(N);
            end
            
            fileList = dir(['N' N '*.mat']);
            fileList = {fileList.name};
        end

        count_done = 0;
        len = length(fileList);
        
        if ~isempty(existsInMat)
            for i = 1:len
                if existInMatfile(fileList{1,i},existInMat)
                    count_done = count_done + 1;
                end
            end

            percent_done = count_done*100/len;
        end
        
        
        if monte
            for i = 1:len
                data = matfile(fileList{1,i});
                [~,s] = size(data.stepInd);
                count_done = count_done + data.stepInd(1,s); 
            end

            percent_done = count_done*100/len*500000;
        end
        
end