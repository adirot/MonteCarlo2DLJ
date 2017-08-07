function unifiedobj = unifyMCobj(objects)
% unify many MC objects to one. objects is a list of MC objects,
% that must have the same simulation parameters. the output object
% is a unified object that has all the steps as if they were all
% preformed in one run.

unifiedobj = objects(1).copyOutputFile('unified');
N = unifiedobj.simulationParam.N;

for i = 1:length(objects)
    objectsi_indIndata = objects(i).indIndata;
    unifiedobj_indIndata = length(unifiedobj.data.allU);
    
    if ~isempty(unifiedobj.data.RDFhisto)
        [~,sizeRDF,~] = size(unifiedobj.data.RDFhisto);
        unifiedobj.data.RDFhisto(1,1:sizeRDF,(unifiedobj_indIndata+1):(unifiedobj_indIndata+1+objectsi_indIndata-3)) =...
            objects(i).data.RDFhisto(1,1:sizeRDF,3:objectsi_indIndata);
    end
    
    if unifiedobj.simulationParam.angleDependent
        unifiedobj.data.allAlphas(1:N,1:N,(unifiedobj_indIndata+1):(unifiedobj_indIndata+1+objectsi_indIndata-3)) =...
            objects(i).data.allAlphas(1:N,1:N,3:objectsi_indIndata);
        unifiedobj.data.allThetas(1:N,1:N,(unifiedobj_indIndata+1):(unifiedobj_indIndata+1+objectsi_indIndata-3)) =...
            objects(i).data.allThetas(1:N,1:N,3:objectsi_indIndata);
        unifiedobj.data.allAngs(1,1:N,(unifiedobj_indIndata+1):(unifiedobj_indIndata+1+objectsi_indIndata-3)) =...
            objects(i).data.allAngs(1,1:N,3:objectsi_indIndata);
    end
    
    unifiedobj.data.allCoords(1:2,1:N,(unifiedobj_indIndata+1):(unifiedobj_indIndata+1+objectsi_indIndata-3)) =...
            objects(i).data.allCoords(1:2,1:N,3:objectsi_indIndata);
    
    try
        unifiedobj.data.allDists(1:N,1:N,(unifiedobj_indIndata+1):(unifiedobj_indIndata+1+objectsi_indIndata-3)) =...
            objects(i).data.allDists(1:N,1:N,3:objectsi_indIndata);
    catch
        [~,~,sizeDists] = size(unifiedobj.data.allDists);
        [~,~,sizeDistsi] = size(objects(i).data.allDists);
        unifiedobj.data.allDists(1:N,1:N,sizeDists) =...
            objects(i).data.allDists(1:N,1:N,sizeDistsi);
    end
    
    if ~isempty(unifiedobj.data.allP)
        unifiedobj.data.allP(1,(unifiedobj_indIndata+1):(unifiedobj_indIndata+1+objectsi_indIndata-3)) =...
            objects(i).data.allP(1,3:objectsi_indIndata);
    end
    
    
    if ~isempty(unifiedobj.data.allPlrc)
        unifiedobj.data.allPlrc(1,(unifiedobj_indIndata+1):(unifiedobj_indIndata+1+objectsi_indIndata-3)) =...
            objects(i).data.allPlrc(1,3:objectsi_indIndata);
    end

    if ~isempty(unifiedobj.data.allV)
        unifiedobj.data.allV(1,(unifiedobj_indIndata+1):(unifiedobj_indIndata+1+objectsi_indIndata-3)) =...
            objects(i).data.allV(1,3:objectsi_indIndata);
    end

    unifiedobj.data.allU(1,(unifiedobj_indIndata+1):(unifiedobj_indIndata+1+objectsi_indIndata-3)) =...
        objects(i).data.allU(1,3:objectsi_indIndata);
    
    if ~isempty(unifiedobj.data.allUlrc)
        unifiedobj.data.allUlrc(1,(unifiedobj_indIndata+1):(unifiedobj_indIndata+1+objectsi_indIndata-3)) =...
            objects(i).data.allUlrc(1,3:objectsi_indIndata);
    end
    
    unifiedobj.data.moveCount =...
        unifiedobj.data.moveCount + objects(i).moveCount;
    
    last_sweepInd = unifiedobj.data.sweepInd(1,unifiedobj_indIndata);
    newsweepInd = [unifiedobj.data.sweepInd(1,1:unifiedobj_indIndata)...
         (objects(i).data.sweepInd(1,3:objectsi_indIndata) + last_sweepInd)];
    unifiedobj.data.sweepInd = newsweepInd;
    
    [~, unifiedobj.data.indIndata] = size(unifiedobj.data.allU); 
    
end

indIndata = unifiedobj.data.indIndata;
unifiedobj.moveCount = unifiedobj.data.moveCount;
unifiedobj.currentCoords = unifiedobj.data.allCoords(1:2,1:N,indIndata);
try
    unifiedobj.currentDists = unifiedobj.data.allDists(1:N,1:N,indIndata);
catch
    [~,~,distsSize] = size(unifiedobj.data.allDists);
    unifiedobj.currentDists = unifiedobj.data.allDists(1:N,1:N,distsSize);
end

unifiedobj.currentU = unifiedobj.data.allU(1,indIndata);

unifiedobj.currentSweep = unifiedobj.data.sweepInd(1,indIndata);

if ~isempty(unifiedobj.data.allP)
    unifiedobj.currentPressure = unifiedobj.data.allP(1,indIndata);
end
   
if ~isempty(unifiedobj.data.allV)
    unifiedobj.currentVir = unifiedobj.data.allV(1,indIndata);
end

if ~isempty(unifiedobj.data.allUlrc)
    unifiedobj.Ulrc = unifiedobj.data.allUlrc(1,indIndata);
end


if ~isempty(unifiedobj.data.allPlrc)
    unifiedobj.Plrc = unifiedobj.data.allPlrc(1,indIndata);
end


if ~isempty(unifiedobj.data.allAlphas)
    unifiedobj.currentAlphas = unifiedobj.data.allAlphas(1:N,1:N,indIndata);
end


if ~isempty(unifiedobj.data.allThetas)
    unifiedobj.currentThetas = unifiedobj.data.allThetas(1:N,1:N,indIndata);
end


if ~isempty(unifiedobj.data.allAngs)
    unifiedobj.currentAngs = unifiedobj.data.allAngs(1,1:N,indIndata);
end

unifiedobj.indIndata = indIndata;
end

