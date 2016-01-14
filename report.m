function percent_done = report(varargin)

        p = inputParser();
        addOptional(p, 'N', ''); 
        addOptional(p, 'fileList', []);
        parse(p, varargin{:});
        Results = p.Results;
        N = Results.N;
        fileList = Results.fileList;
        
        if isempty(fileList)
            if ~strcmp(N,'')
                N = num2str(N);
            end
            
            fileList = dir(['N' N '*.mat']);
            fileList = {fileList.name};
        end

        count_done = 0;
        len = length(fileList);
        
        for i = 1:len
            if existInMatfile(fileList{1,i},'RDFhisto')
                count_done = count_done + 1;
            end
        end
        
        percent_done = count_done*100/len;
end