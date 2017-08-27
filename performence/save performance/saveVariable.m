function [timeOfEachStep, steps] = saveVariable(var,fileName,varargin)
% performence cheack for matlab's save function.


    p = inputParser();
    % matfile version
    addOptional(p, 'version','-7.3');
    % if the variable size in the matfile should be updaetd each time we
    % save: true. if the the varible size in the matfile should be set
    % upfront: false. (true can be used only for v7.3)
    addOptional(p, 'updateVarSize', false);
    % how many steps to save each time we save. saveBlocSize should be a
    % devisor of the number of steps steps
    addOptional(p, 'saveBlocSize', 1);
    parse(p,varargin{:});
    results = p.Results;
    version = results.version;
    updateVarSize = results.updateVarSize;
    saveBlocSize = results.saveBlocSize;
    
    sizeVar = size(var);
    N = sizeVar(2);
    numOfSteps = sizeVar(3);
    
    % Check input 
    if mod(numOfSteps,saveBlocSize) > 0 
        error('saveBlock size should be a devisor of the number of steps. (the size of var is 1 by N by steps)');
    end
    
    switch version
        
        % version 7.3 allows to load and save parts of a variable, so we
        % use that ability.
        case '-v7.3' 
            matFilev7_3 = matfile(fileName);
            matFilev7_3 = matfile(fileName,'Writable',true);
            
            
                
            if updateVarSize
                s = zeros(2,N,2);
                matFilev7_3.var = s;
                steps = numOfSteps;
            else
                s = zeros(sizeVar);
                matFilev7_3.var = s;
            end
            counter = 0;
            for i = 1:saveBlocSize:numOfSteps
                tic;
                matFilev7_3.var(:,:,i:(i+saveBlocSize-1))...
                    = var(:,:,i:(i+saveBlocSize-1));
                tocVal = toc;
                counter = counter + 1;
                timeOfEachStep(1,counter) = tocVal;
                steps(1,counter) = (i+saveBlocSize-1);
                disp(i);
            end
        
        % version 7 does not allow to save an load parts of a variable, so
        % we use append to create another variable for each bloc we save
        case {'-v7', '-v6'}
            counter = 0;
            for i = 1:saveBlocSize:numOfSteps
                tic;
                varBlocName = ['bloc' num2str(i) 'to'...
                    num2str(i+saveBlocSize-1)];
                assignin('base', varBlocName,...
                    var(:,:,i:(i+saveBlocSize-1)));
                save(fileName,varBlocName,'-append',version);
                counter = counter + 1;
                tocVal = toc;
                timeOfEachStep(1,counter) = tocVal;
                steps(1,counter) = (i+saveBlocSize-1);
                disp(i);
            end   
    end

end