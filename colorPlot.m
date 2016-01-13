function colorPlot(x,y,varargin)
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
addOptional(p, 'lineStyle', '+');
addOptional(p, 'colormap','jet');
addOptional(p, 'figHandle', []);


parse(p, varargin{:});
Results = p.Results;

addLegend = Results.addLegend;
plotBar = Results.plotBar;
Interpreter = Results.Interpreter;
length2plot = Results.length2plot;
style = Results.lineStyle;
colormap = Results.colormap;
figHandle = Results.figHandle;


    [numOfPlots,~] = size(x);
    
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
       figHandle = figure;
   else
       figure(figHandle);
   end
        
   hold on;   
   for i = 1:numOfPlots
       if plotBar
            increment = x(i,2)-x(i,1);
            if isempty(length2plot)
                bar(x(i,:)+increment/2,y(i,:),1,'EdgeColor',cmap(i,:)...
                    ,'FaceColor','none')
            else
                if length2plot(i) > 0
                    bar(x(i,1:length2plot(i))+increment/2,...
                        y(i,1:length2plot(i)),1,'EdgeColor',cmap(i,:)...
                        ,'FaceColor','none')
                end
            end
       else
           if isempty(length2plot)
               
                    plot(x(i,:),y(i,:),'color',cmap(i,:),'lineStyle',style);
               

           else
                if length2plot > 0
                    
                    plot(x(i,1:length2plot(i)),y(i,1:length2plot(i)),...
                        'color',cmap(i,:),'lineStyle',style);
                    
                    
                end
           end
               
       end
   end
   
   if ~isempty(addLegend)
        legend(addLegend,'Interpreter', Interpreter)
   end
   
end
       