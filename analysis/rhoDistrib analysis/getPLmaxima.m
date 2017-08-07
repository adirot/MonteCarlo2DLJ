% find PL maxima
N = 625;
sq = 256;
m1 = 3;
algo = '';
% 'b':0.1 'r':0.2 'k':0.01 'g':0.05
col = 'b';
rhostr = '0_5';
rho = num2str(strrep(rhostr,'_','.'));
%list = dir(['r*' algo '*rho' rhostr '*m' num2str(m1) '*' num2str(sq) '.fig']);
list = dir('*T0_0*rho0_2*m3*121sys.fig');
list = {list.name};
clear maxima;

for i = 1:length(list)
    name = list{1,i};
    locT = strfind(name,'T');
    %locrho = strfind(name,'rho0');
    %locm = strfind(name,'m');
    %locm = locm(1);
    T = str2double(strrep(name((locT+1):(locrho-1)),'_','.'));
    rho = 0.2; %rho = str2double(strrep(name((locrho+4):(locm-1)),'_','.'));
    m = 3; %m = str2double(name(locm+1)); 
    snap = dir(['s*T' strrep(num2str(T),'.','_') 'rho' strrep(num2str(rho),'.','_') '*m' num2str(m) '*.fig']);
    snap = {snap.name}; 
%    if length(snap) ~= 1
%         if length(snap) > 1
%             error('snap > 1');
%         end
%    else
        try
            f = open(snap{1,1});
            movegui(f,[100 100]);
        catch
        end
 %   end
    
    open(list{1,i}); hold on; xlim([0 10]);
    h = gcf; %current figure handle
    axesObjs = get(h, 'Children');  %axes handles
    dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
    objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
    xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
    ydata = get(dataObjs, 'YData');
    [pks,locs] = findpeaks(ydata);
        plot(xdata(locs),pks,'+r');
        [xuser,yuser] = ginput(2);
        maxima{i,1} = xuser;
        maxima{i,2} = yuser;
    
    close all;
end

xmaxima = []; 
ymaxima = [];
for i = 1:length(maxima)
    if ~isempty(maxima{i,1})
        numOfMax = length(maxima{i,1});
        xmaxima = [xmaxima maxima{i,1}'];
        name = list{1,i};
        locT = strfind(name,'T');
        locrho = strfind(name,'rho0');
        locm = strfind(name,'m');
        locm = locm(1);
        T = str2double(strrep(name((locT+1):(locrho-1)),'_','.'));
        ymaxima = [ymaxima T];
        if numOfMax == 2
            ymaxima = [ymaxima T];
        end
        disp(i);
    end
end
figure;
plot(xmaxima*rho,ymaxima,'color',col,'Marker','+','lineStyle','none');
%herrorbar(xmaxima*rho,ymaxima,sq*ones(1,length(xmaxima))*rho/(N*2),'+');
hold on;

% mean field :
f = @(y) (y.*cosh(y)-sinh(y))./(sinh(y).*cosh(y)-y);
c = @(s) (cosh(s));
g = @(s) (1 + 2.*c(s).*f(s) + (f(s)).^2);

T = @(s) 27*f(s/2).*(f(s/2)+c(s/2))./(4*(g(s/2)).^2); % Normalized by Tc = 8a/(27b)
P = @(s) 27*((f(s/2)).^2).*(1-(f(s/2)).^2)./((g(s/2)).^2); % Normalized by Pc = a / (27b^2)
Vmean = @(s) 3*f(s/2).*(c(s/2)+f(s/2))./g(s/2); % Normalized by Vc = 3b
deltaV = @(s) 6*f(s/2).*sinh(s/2)./g(s/2); % Normalized by Vc

s = 0:0.1:100;
m=m1;

% scr:
% plot((2/(3*pi))*((12/m)^(2/(m-12)))*(Vmean(s)+deltaV(s)/2),(16/27)*((4/(m-2))*((12/m)^(m/(m-12)))-(4/10)*((12/m)^(12/(m-12))))*T(s),'k');
% plot((2/(3*pi))*((12/m)^(2/(m-12)))*(Vmean(s)-deltaV(s)/2),(16/27)*((4/(m-2))*((12/m)^(m/(m-12)))-(4/10)*((12/m)^(12/(m-12))))*T(s),'k');

% hcr:
plot((Vmean(s)-deltaV(s)/2)/(6*pi),(2^(2-m))*(4/(27*(m-2)))*T(s),'b');
plot((Vmean(s)+deltaV(s)/2)/(6*pi),(2^(2-m))*(4/(27*(m-2)))*T(s),'b');

% E Vs rho
% figure;
% alpha = -1;
% m = 3;
% B = 1;
% C = 1;
% t = 1;

% scr:
% E_op6 = @(B,C,t,m,alpha,s) ((16/27)*((C/(m-2))*((12*B/(m*C))^(m/(m-12)))-(B/10)*(12*B/(m*C))^(12/(m-12)))*(T(s)/t)).^((m-12)/(12*alpha)); 
% rho_op6 = @(B,C,t,m,alpha,s) ((2/(3*pi))*((12*B/(m*C))^(2/(m-12)))*((E_op6(B,C,t,m,alpha,s)).^(-2*alpha/(m-12))));
% plot(rho_op6(B,C,t,m,alpha,s).*(Vmean(s)+deltaV(s)/2),E_op6(B,C,t,m,alpha,s));
% hold on;
% plot(rho_op6(B,C,t,m,alpha,s).*(Vmean(s)-deltaV(s)/2),E_op6(B,C,t,m,alpha,s));
% 
% E_op6_sim = @(B,C,t,m,alpha,ymaxima) ((B/C)^(1/alpha))*((4*t/B)^((12-m)/(12*alpha)))*(ymaxima.^((m-12)/(12*alpha)));
% rho_op6_sim = @(B,C,t,m,alpha,ymaxima) ((C/B)^(2/(12-m)))*(E_op6_sim(B,C,t,m,alpha,ymaxima).^(2*alpha/(12-m)));
% plot(rho_op6_sim(B,C,t,m,alpha,Trhom3sq256(2,:)).*Trhom3sq256(1,:),E_op6_sim(B,C,t,m,alpha,Trhom3sq256(2,:)));

% hcr:
% rc = 1;
% t = 1;
% B = 1;
% C = 1;
% m = 6;
% alpha = 1;
% figure;
% E_op6_hcr = @(rc,C,t,m,alpha,s) (27*(m-2)*(rc^m)*t./(4*C*T(s)*(2^(2-m)))).^(1/alpha); 
% rho_op6_hcr = @(rc) (1/(6*pi*(rc^2)));
% plot(rho_op6_hcr(rc)*(Vmean(s)+deltaV(s)/2),E_op6_hcr(rc,C,t,m,alpha,s));
% % normalized:
% %plot((Vmean(s)+deltaV(s)/2),E_op6_hcr(rc,C,t,m,alpha,s)/E_op6_hcr(rc,C,t,m,alpha,0.0001));
% hold on;
% plot(rho_op6_hcr(rc)*(Vmean(s)-deltaV(s)/2),E_op6_hcr(rc,C,t,m,alpha,s));
% %plot((Vmean(s)-deltaV(s)/2),E_op6_hcr(rc,C,t,m,alpha,s)/E_op6_hcr(rc,C,t,m,alpha,0.0001));
% 
% % E_op6_sim_hcr = @(rc,C,t,m,alpha,ymaxima) (t*(rc^m)./(C*ymaxima)).^(1/alpha);
% rho_op6_sim_hcr = @(rc) 1/rc^2;
% %plot(rho_op6_sim_hcr(rc).*m3Psys256hcr(:,1),E_op6_sim_hcr(rc,C,t,m,alpha,m3Psys256hcr(:,2)));
% a1 = 12; b1 = 10; 
% plot(a1*a(:,1),b1*E_op6_sim_hcr(rc,C,t,m,alpha,a(:,2)),'+');
