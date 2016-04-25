function [fitresult, mfit, mError, gof] =...
    createFitRDF(xs, ys, T, rho,m , steps, varargin)
%CREATEFIT(X,Y)
%  Create a fit for RDF.
%
%  Data for 'untitled fit 1' fit:
%      X Input : x
%      Y Output: y
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

p = inputParser();
addOptional(p, 'plotFig', true);
parse(p, varargin{:});
Results = p.Results;
plotFig = Results.plotFig;
                    

[Nplots, ~] = size(ys);

if plotFig
    h = colorPlot(xs,ys,'addLegend',steps);
    hold on;
    title(['RDF with fit for T = ' num2str(T) ' \rho = '...
        num2str(rho) ' m = ' num2str(m)]);
    message = sprintf(['4*(1/' num2str(T) ')*((1/x)^{12} - (1/x)^m)']);
    text(2,7,message);
end

% Set up fittype and options.
ft = fittype( ['4*(1/' num2str(T) ')*((1/x)^12 - (1/x)^m)'],...
        'independent', 'x', 'dependent', 'y' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = -Inf;
opts.StartPoint = 0.970592781760616;
opts.Upper = Inf;


for i = 1:Nplots
    y = ys(i,:);
    x = xs(i,:);
    [xData, yData] = prepareCurveData( x, y );


    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    conf = confint(fitresult);
    mfit(i) = coeffvalues(fitresult);
    mError(i,1) = conf(1);
    mError(i,2) = conf(2);
    
    if plotFig
        % Plot fit 
        plot( fitresult );
        %legend( h, 'simulation', 'fit', 'Location', 'NorthEast' );
        
        % annotation(gcf,'textbox',...
        %     [0.254571428571429 0.654761904761905 0.181142857142857 0.20952380952381],...
        %     'String',{[formula(fitresult) '\n' 'm = ' coeffvalues(fitresult)...
        %     ' (' num2str(conf(1)) ',' num2str(conf(2)) ')']});
        message = sprintf(['steps: ' steps{i} '\n' 'm = '...
            num2str(coeffvalues(fitresult))...
            ' (' num2str(conf(1)) ',' num2str(conf(2)) ')']);
        text(2,7-i,message,'EdgeColor','k');
    end
end

if plotFig
    ylim([-3 8])
    set(findall(gcf,'-property','FontSize'),'FontSize',20);

    % Label axes
    xlabel( 'distance' );
    ylabel( '-log(RDF)' );
end
