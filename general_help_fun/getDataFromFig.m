function [xdata,ydata] = getDataFromFig(fileName)

% Extract data from figure (only supports line data)
openfig(fileName,'new','invisible');
h = gcf;
axesObjs = get(h,'Children');
dataObjs = get(axesObjs,'Children');
%objTypes = get(dataObjs,'Type');
lineData = findobj(dataObjs,'Type','line');
xdata = get(lineData,'XData');
ydata = get(dataObjs,'YData');


end

