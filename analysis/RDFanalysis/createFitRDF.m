function [fitresult, mfit, mError, nfit, nError, Tfit, TError, gof] =...
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
addOptional(p, 'plotFig', false);
addOptional(p, 'freen', false);
addOptional(p, 'freeTandn', false);
addOptional(p, 'freeTnbound', false);
addOptional(p, 'freeTnset', false);
addOptional(p, 'nset', []);
addOptional(p, 'hcr', false);
addOptional(p, 'hcr_set_T', false);
addOptional(p, 'ignoreZerosAtend', false);
parse(p, varargin{:});
Results = p.Results;
plotFig = Results.plotFig;
freen = Results.freen;
freeTandn = Results.freeTandn;
freeTnbound = Results.freeTnbound;
freeTnset = Results.freeTnset;
nset = Results.nset; 
hcr = Results.hcr;
hcr_set_T = Results.hcr_set_T;
ignoreZerosAtend = Results.ignoreZerosAtend;

fitresult=[]; mfit=[]; mError=[]; nfit=[]; nError=[]; Tfit=[]; TError=[]; gof=[];
fitProblem = false;

if freeTandn
    freen = false;
end

if and(freeTnset,isempty(nset))
    error('choose n with nset input parameter');
end

[Nplots, ~] = size(ys);

if plotFig
    h = colorPlot(xs,ys,'addLegend',steps);
    hold on;
    if freen
        title(['RDF with fit for T = ' num2str(T) ' \rho = '...
            num2str(rho) ' m = ' num2str(m) ' n free']);
        message = sprintf(['4*(1/' num2str(T) ')*((1/x)^n - (1/x)^m)']);
    else
        if freeTandn
            title(['RDF with fit for T = ' num2str(T) ' \rho = '...
                num2str(rho) ' m = ' num2str(m) ' n,T free']);

            message = sprintf('4*(1/T)*((1/x)^n - (1/x)^m)');
        
        else
            if freeTnbound
                title(['RDF with fit for T = ' num2str(T) ' \rho = '...
                    num2str(rho) ' m = ' num2str(m) ' n set bound, T free']);

                message = sprintf(['4*(1/' num2str(T) ')*((1/x)^{12} - (1/x)^m)']);
            
            else
                if freeTnset

                    title(['RDF with fit for T = ' num2str(T) ' \rho = '...
                        num2str(rho) ' m = ' num2str(m) ' n set to ' num2str(nset)]);

                    message = sprintf(['4*(1/T)*((1/x)^{' num2str(nset) '} - (1/x)^m)']); 
                else
                    title(['RDF with fit for T = ' num2str(T) ' \rho = '...
                        num2str(rho) ' m = ' num2str(m) ' n set to 12']);

                    message = sprintf(['4*(1/' num2str(T) ')*((1/x)^{12} - (1/x)^m)']);
                    if hcr
                        title(['RDF with fit for T = ' num2str(T) ' \rho = '...
                        num2str(rho) ' m = ' num2str(m) ' hcr ']);

                        message = sprintf('-(1/T)*((x)^(-m))');
                    else
                        if hcr_set_T
                            title(['RDF with fit for T = ' num2str(T) ' \rho = '...
                            num2str(rho) ' m = ' num2str(m) ' hcr with set T']);

                            message = sprintf(['(1/' num2str(T) ')*(- (1/x)^m)']);
                        end 
                    end
                end
            end
        end
    end
    
    text(2,7,message);
end

% Set up fittype and options.
if freen
    ft = fittype( ['4*(1/' num2str(T) ')*((1/x)^n - (1/x)^m)'],...
        'independent', 'x', 'dependent', 'y' );
else
    if freeTandn
        ft = fittype('4*(1/T)*((1/x)^n - (1/x)^m)',...
                'independent', 'x', 'dependent', 'y' );
    
    else
        if freeTnbound
            ft = fittype( '4*(1/T)*((1/x)^n - (1/x)^m)',...
                    'independent', 'x', 'dependent', 'y' );
            opts = fitoptions( ft );
            opts.Display = 'Off';
            
            % set the bounds of n to 10-20
            opts.Lower = [-Inf -Inf 10]; 
            opts.Upper = [Inf Inf 20];
        else
            if freeTnset
                ft = fittype( ['4*(1/T)*((1/x)^' num2str(nset) ' - (1/x)^m)'],...
                        'independent', 'x', 'dependent', 'y' );
            else
                ft = fittype( ['4*(1/' num2str(T) ')*((1/x)^12 - (1/x)^m)'],...
                        'independent', 'x', 'dependent', 'y' );
                if hcr
                    ft = fittype('-(1/T)*((x)^(-m))',...
                        'independent', 'x', 'dependent', 'y' );
                else
                    if hcr_set_T
                        ft = fittype(['(1/' num2str(T) ')*(- (1/x)^m)'],...
                            'independent', 'x', 'dependent', 'y' );
                    end 
                end
            end
        end
    end
end
% opts = fitoptions( ft );
% opts.Display = 'Off';
% opts.Lower = -Inf;
% opts.StartPoint = 0.970592781760616;
% opts.Upper = Inf;

nfit = []; nError = [];
Tfit = []; TError = [];
mfit = []; mError = [];

for i = 1:Nplots
    
    %  hcr should only be done for r > 2. for 300 bins and rcutoff 10 its
    %  from 32
    if hcr
        [~, hcrInd] = min(abs(xs(i,:)-2));
        sind = hcrInd + 2;
    else
        sind = 1;
    end
    
    if ignoreZerosAtend
        lind = max(find(xs(i,:)));
        lind = lind - 1;
    else
        lind = 299;
    end
    
    y = ys(i,sind:lind);
    x = xs(i,sind:lind);
    [xData, yData] = prepareCurveData( x, y );

    % stop warning about no starting point provided
    warning('off','curvefit:fit:noStartPoint');
    % stop warning about removing nan and inf data
    warning('off','curvefit:prepareFittingData:removingNaNAndInf');

    % Fit model to data.
%     [fitresult, gof] = fit( xData, yData, ft, opts );
     try
        [fitresult, gof] = fit( xData, yData, ft);
        conf = confint(fitresult);
        coeff = coeffvalues(fitresult);

        mfit{i} = coeff(1);
        mError{i,1} = conf(1,1);
        mError{i,2} = conf(2,1);
        if freen 
            nfit{i} = coeff(2);
            nError{i,1} = conf(1,2);
            nError{i,2} = conf(2,2);
        else
            if or(freeTandn,freeTnbound)
                Tfit{i} = coeff(1);
                TError{i,1} = conf(1,1);
                TError{i,2} = conf(2,1);
                mfit{i} = coeff(2);
                mError{i,1} = conf(1,2);
                mError{i,2} = conf(2,2);
                nfit{i} = coeff(3);
                nError{i,1} = conf(1,3);
                nError{i,2} = conf(2,3);
            end

            if freeTnset
                Tfit{i} = coeff(1);
                TError{i,1} = conf(1,1);
                TError{i,2} = conf(2,1);
                mfit{i} = coeff(2);
                mError{i,1} = conf(1,2);
                mError{i,2} = conf(2,2);
            end
        
                
        end
        
        if hcr
                Tfit{i} = coeff(1);
                TError{i,1} = conf(1,1);
                TError{i,2} = conf(2,1);
                mfit{i} = coeff(2);
                mError{i,1} = conf(1,2);
                mError{i,2} = conf(2,2);
        end
        
        if hcr_set_T
                Tfit{i} = coeff(1);
                TError{i,1} = conf(1,1);
                TError{i,2} = conf(2,1);
                mfit{i} = coeff(2);
                mError{i,1} = conf(1,2);
                mError{i,2} = conf(2,2);
        end
            
    catch ME
        conf = [];
        coeff = [];

        mfit{i} = [];
        mError{i,1} = [];
        mError{i,2} = [];
        
        if freen 
            nfit{i} = [];
            nError{i,1} = [];
            nError{i,2} = [];
        else
            if or(freeTandn,freeTnbound)
                Tfit{i} = [];
                TError{i,1} = [];
                TError{i,2} = [];
                mfit{i} = [];
                mError{i,1} = [];
                mError{i,2} = [];
                nfit{i} = [];
                nError{i,1} = [];
                nError{i,2} = [];
            end

            if freeTnset
                Tfit{i} = [];
                TError{i,1} = [];
                TError{i,2} = [];
                mfit{i} = [];
                mError{i,1} = [];
                mError{i,2} = [];
            end
        
            if hcr
                Tfit{i} = [];
                TError{i,1} = [];
                TError{i,2} = [];
                mfit{i} = [];
                mError{i,1} = [];
                mError{i,2} = [];
            end
            if hcr_set_T
                Tfit{i} = [];
                TError{i,1} = [];
                TError{i,2} = [];
                mfit{i} = [];
                mError{i,1} = [];
                mError{i,2} = [];
            end
                
        end        
            disp('fit problem!');
            fitProblem = true;
            disp(ME.message);
    end
end
    
    if and(plotFig,~fitProblem)
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
            if or(freeTandn,freeTnbound)
                message = sprintf(['steps: ' steps{i} '\n' 'T = '...
                    num2str(coeff(1))...
                    ' (' num2str(conf(1,1)) ',' num2str(conf(2,1)) ')\nn = '...
                    num2str(coeff(2))...
                    ' (' num2str(conf(1,2)) ',' num2str(conf(2,2)) ')\nm = '...
                    num2str(coeff(3))...
                    ' (' num2str(conf(1,3)) ',' num2str(conf(2,3)) ')']);
            else
                if freeTnset
                    message = sprintf(['steps: ' steps{i} '\n' 'T = '...
                        num2str(coeff(1))...
                        ' (' num2str(conf(1,1)) ',' num2str(conf(2,1)) ')\nm = '...
                        num2str(coeff(2))]);
                else
                    
                    message = sprintf(['steps: ' steps{i} '\n' 'm = '...
                        num2str(coeff(1))...
                        ' (' num2str(conf(1,1)) ',' num2str(conf(2,1)) ')']);
                end
            end
        end
        text(2,7-i,message,'EdgeColor','k');
    end
end

% if plotFig
%     ylim([-3 8])
%     set(findall(gcf,'-property','FontSize'),'FontSize',12);
% 
%     % Label axes
%     xlabel( 'distance' );
%     ylabel( '-log(RDF)' );
% end
