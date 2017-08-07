function [xdata,ydata] = getDataFromFig(fileName)

% Extract data from figure (only supports line data)
openfig(fileName,'new','invisible');
h = gcf;
axesObjs = get(h,'Children');
dataObjs = get(axesObjs,'Children');
%objTypes = get(dataObjs,'Type');
try 
    lineData = findobj(dataObjs,'Type','line');
catch
        lineData = findobj(dataObjs{2,1},'Type','line');
end
xdata = get(lineData,'XData');
ydata = get(lineData,'YData');


end

