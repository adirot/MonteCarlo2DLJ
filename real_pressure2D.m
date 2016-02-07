function [P,rho,rho_times_B2,good_ind,cantusevirial_ind...
    ,problem_ind] = real_pressure2D(T,rho,my_eps,potantial,varargin)

%% find the pressure of a Lennard-Jonse 2 dimantional gas with reduced
%% temprature T and reduced density rho, up to the second virial term.

% The virial expantion for a gas of 2 dimantional interacting particals:
% (in reduced units)

%       P = T*rho(1 - rho*B2)

%       where

%       B2 = pi * int(exp(-4 beta((1/r)^12 - (1/r)^6)) - 1)rdr

% The analitical solution is:

%       B2 = pi*(1/(2^(2/3)))beta^(1/6) (-Gamma[5/6] F1[-(1/6), 1/2, beta] + 
%               Sqrt[beta] Gamma[4/3] F1[1/3, 3/2, beta])

% where gamma is the gamma function, and F1 is the 
%   Kummer confluent hypergeometric function  
% (http://mathworld.wolfram.com/ConfluentHypergeometricFunctionoftheFirstKind.html)

% some values calculated with wolfram alpha:
% Tdb = [0.1 0.2 0.3 0.4 0.45 0.5 0.6 0.7 0.9 1 10 100];
% B2db = [8948.07, 95.165, 22.1906, 10.192, 7.65916, 5.98671, 3.94465, ...
%     2.76122, 1.46388, 1.07347, -1.08036, -0.951233];
% B2stdb = [1365.01, 21.4645 , 5.40892, 2.25374 , 1.52603, 1.0259,...
%     0.387399, -0.00116176, -0.448385 , -0.588563 , -1.49065 , -1.56293];

p = inputParser();
addOptional(p,'m',6);
parse(p, varargin{:});
Results = p.Results;
m = Results.m;

P = zeros(1,length(rho));
rho_times_B2 = zeros(1,length(rho));
good_ind = [];
cantusevirial_ind = [];
problem_ind = [];

for ii = 1:length(rho)
    [P(ii),errorind,rho_times_B2(ii)] =...
        calc_real_pressure(T,rho(ii),my_eps,potantial,m);
    
    switch errorind
        case 0 % no error
            good_ind = [good_ind ii];
        case 1 % |virial| > 1
            cantusevirial_ind = [cantusevirial_ind ii];
        case 2 % 0.5 < |virial| < 1
            problem_ind  = [problem_ind ii];
    end
            
end

    function [P,errorind,rho_times_B2] = ...
            calc_real_pressure(T,rho,my_eps,potantial,m)
        
        switch potantial
            case 'st'
                B2 = B2stdb(ind);
            otherwise
                B2 = calc_B2(T,m);
        end


            rho_times_B2 = rho*B2;
            second_term = T*rho*rho_times_B2;
            first_term = T*rho;

                if abs(first_term) < abs(second_term)   
                    P = 0;
                    errorind = 2;
                else 
                    P = first_term - second_term;

                    if abs(second_term)/abs(first_term) < my_eps
                        errorind = 0;
                    else 
                        errorind = 1;
                    end
                end
                
    end
    

 
end

