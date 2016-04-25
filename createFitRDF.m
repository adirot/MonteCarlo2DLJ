function [fitresult, mfit, mError, nfit, nError, gof] =...
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
addOptional(p, 'freen', false);
parse(p, varargin{:});
Results = p.Results;
plotFig = Results.plotFig;
freen = Results.freen;
                    

[Nplots, ~] = size(ys);

if plotFig
    h = colorPlot(xs,ys,'addLegend',steps);
    hold on;
    if freen
        title(['RDF with fit for T = ' num2str(T) ' \rho = '...
            num2str(rho) ' m = ' num2str(m) ' n free']);
        message = sprintf(['4*(1/' num2str(T) ')*((1/x)^n - (1/x)^m)']);
    else
        title(['RDF with fit for T = ' num2str(T) ' \rho = '...
            num2str(rho) ' m = ' num2str(m) ' n set to 12']);
        
        message = sprintf(['4*(1/' num2str(T) ')*((1/x)^{12} - (1/x)^m)']);
    end
    
    text(2,7,message);
end

% Set up fittype and options.
if freen
    ft = fittype( ['4*(1/' num2str(T) ')*((1/x)^n - (1/x)^m)'],...
        'independent', 'x', 'dependent', 'y' );
else
    ft = fittype( ['4*(1/' num2str(T) ')*((1/x)^12 - (1/x)^m)'],...
            'independent', 'x', 'dependent', 'y' );
end
% opts = fitoptions( ft );
% opts.Display = 'Off';
% opts.Lower = -Inf;
% opts.StartPoint = 0.970592781760616;
% opts.Upper = Inf;


for i = 1:Nplots
    y = ys(i,:);
    x = xs(i,:);
    [xData, yData] = prepareCurveData( x, y );


    % Fit model to data.
%     [fitresult, gof] = fit( xData, yData, ft, opts );
    [fitresult, gof] = fit( xData, yData, ft);
    conf = confint(fitresult);
    coeff = coeffvalues(fitresult);
    
    mfit(i) = coeff(1);
    mError(i,1) = conf(1,1);
    mError(i,2) = conf(2,1);
    if freen
        nfit(i) = coeff(2);
        nError(i,1) = conf(1,2);
        nError(i,2) = conf(2,2);
    else 
        nfit = []; nError = [];
    end
        
    if plotFig
        % Plot fit 
        plot( fitresult );
        %legend( h, 'simulation', 'fit', 'Location', 'NorthEast' );
        
        % annotation(gcf,'textbox',...
        %     [0.254571428571429 0.654761904761905 0.181142857142857 0.20952380952381],...
        %     'String',{[formula(fitresult) '\n' 'm = ' coeffvalues(fitresult)...
        %     ' (' num2str(conf(1)) ',' num2str(conf(2)) ')']});
        if freen
            message = sprintf(['steps: ' steps{i} '\n' 'm = '...
                num2str(coeff(1))...
                ' (' num2str(conf(1,1)) ',' num2str(conf(2,1)) ')\nn = '...
                num2str(coeff(2))...
                ' (' num2str(conf(1,2)) ',' num2str(conf(2,2)) ')']);
        else
            message = sprintf(['steps: ' steps{i} '\n' 'm = '...
                num2str(coeff(1))...
                ' (' num2str(conf(1,1)) ',' num2str(conf(2,1)) ')']);
        end
        text(2,7-i,message,'EdgeColor','k');
    end
end

if plotFig
    ylim([-3 8])
    set(findall(gcf,'-property','FontSize'),'FontSize',16);

    % Label axes
    xlabel( 'distance' );
    ylabel( '-log(RDF)' );
end
