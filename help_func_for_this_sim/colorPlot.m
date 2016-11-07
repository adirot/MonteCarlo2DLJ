function figHandle = colorPlot(x,y,varargin)
% every row in x is an x axis to plot
% every rew in y is a y axis to plot
% this will plot all [x, y]'s in different colors
   
% optional: 'labels': a cell array of strings to be used as legend
% optional: 'plotBar': plot as bar
% otional: 'Interpreter' for labels

p = inputParser();

addOptional(p, 'addLegend', []);
addOptional(p, 'plotBar', false);
addOptional(p, 'Interpreter', 'none');
addOptional(p, 'length2plot', []);
addOptional(p, 'plotFrom', []);
addOptional(p, 'colormap','jet');
addOptional(p, 'figHandle', []);
addOptional(p, 'lineStyle', '-');

parse(p, varargin{:});
Results = p.Results;

addLegend = Results.addLegend;
plotBar = Results.plotBar;
Interpreter = Results.Interpreter;
length2plot = Results.length2plot;
plotFrom = Results.plotFrom;
colormap = Results.colormap;
figHandle = Results.figHandle;
lineStyle = Results.lineStyle;

[numOfPlots,Nx] = size(x);

if isempty(plotFrom)
    plotFrom = ones(1,numOfPlots);
end

if isempty(length2plot)
    length2plot = Nx*ones(1,numOfPlots);
end


    
  switch colormap
      case 'hsv'
            cmap = hsv(numOfPlots);
      
      case 'winter'
          cmap = winter(numOfPlots);
      case 'jet'
          cmap = jet(numOfPlots);
      otherwise
          cmap = jet(numOfPlots);
  end
          
   if isempty(figHandle)
       figHandle = figure('units','normalized','outerposition',[0 0 1 1]);
   else
       figure(figHandle);
   end
        
   hold on;
   
   for i = 1:numOfPlots
       if plotBar
            increment = x(i,2)-x(i,1);
            
            bar(x(i,plotFrom(i):length2plot(i))+increment/2,...
                        y(i,plotFrom(i):length2plot(i)),...
                        1,'EdgeColor',cmap(i,:)...
                        ,'FaceColor','none')
            
       else
           plot(x(i,plotFrom(i):length2plot(i)),...
               y(i,plotFrom(i):length2plot(i)),...
                        'color',cmap(i,:),'lineStyle',lineStyle);
               
       end
   end
   
   if ~isempty(addLegend)
        legend(addLegend,'Interpreter', Interpreter)
   end
   
end
       
