%--------------------------------------------------------------------------
%                          Specify data inputs
%--------------------------------------------------------------------------
close all
clear
fileLoc='C:\....';
fileName='Filename....xlsx';
Delta=readmatrix([fileLoc,fileName],'Sheet','Isotopes','Range','##:##');   % Soil and Plant Water Isotopes
Sap_H2_meas=readmatrix([fileLoc,fileName],'Sheet','Isotopes','Range','##:##');   % Tree Sap Water Delta-2-H
Sap_O18_meas=readmatrix([fileLoc,fileName],'Sheet','Isotopes','Range','##:##');   % Tree Sap Water Delta-18-O
SoilPsi=readmatrix([fileLoc,fileName],'Sheet','Soil_Potential','Range','##:##');   % Soil water potential
SapFlux=readmatrix([fileLoc,fileName],'Sheet','Sap_Flow','Range','##:##');   % Tree sap flux
Roots=readmatrix([fileLoc,fileName],'Sheet','Roots','Range','##:##');   % Root distribution
numVars=3; varNames={'Time','Psi'}; varTypes={'double','double'}; 
dataRange='##:##'; dataSheet='Leaf_Potential'; % Leaf Water Potential
opts = spreadsheetImportOptions('NumVariables',numVars,'VariableNames',varNames, ...
    'VariableTypes',varTypes,'Sheet',dataSheet,'DataRange',dataRange);
LeafPsi=readmatrix([fileLoc,fileName],opts);
%--------------------------------------------------------------------------
LeafPsi(:,1)=LeafPsi(:,1)*24;   % convert units from days to hours.
SapFlux(:,1)=SapFlux(:,1)*24;   % convert units from days to hours.

mm=2500; % Number of intervals to interpolate soil water potential depths and leaf water potential times.
Zr=250; % Maximum root depth.
Xylem_count = ####; % Monte Carlo simulation length. i.e., number of samples to be taken from the distribution of xylem water isotopic compositions (Xylem_dist).
cond = #.####; % Conductance of the flow path, K(t).

leafStime=min(LeafPsi(:,1));    % Start time of leaf potential measurement.
leafEndtime=max(LeafPsi(:,1));  % End time of leaf potential measurement.
psidepth=zeros(mm,1);           % Vector for storing soil water potential depth intervals.
leafPsitime=zeros(mm,1);        % Vector for storing times of leaf water potentials.

for i=1:mm
    psidepth(i)=(i-1)*Zr/(mm-1);
    leafPsitime(i)=leafStime+(i-1)*(leafEndtime-leafStime)/(mm-1);
end

%--------------------------------------------------------------------------
% The following are for interpolating leaf water potential, soil potential, and soil isotope values for integration
%--------------------------------------------------------------------------
Interpolation="matlab";
InterpMethod=["nearestinterp" "pchip" "smoothingspline" "cubicinterp"];
fitType=InterpMethod(1);
% Make sure leaf potential and sap flux measurements have similar start and end time.
% Also make sure soil isotope and soil potenital measurements start at soil surface and end at depth Zr.

C2H_M=csape(Delta(:,1)/Zr,[0,(Delta(:,2)'+1000)/6.420135,0],[1 1]);
C18O_M=csape(Delta(:,1)/Zr,[0,(Delta(:,3)'+1000)*2.0052,0],[1 1]);
soilPsi=csape(SoilPsi(:,1)/Zr,[0,SoilPsi(:,2)',0],[2 2]);
sapFlux=csape(SapFlux(:,1),[0,SapFlux(:,2)', 0],[1 1]);
leafPsi=csape(LeafPsi(:,1),[0, LeafPsi(:,2)', 0],[1 1]);
roots=csape(Roots(:,1)/Zr,[0,Roots(:,2)',0],[2 2]);

%--------------------------------------------------------------------------
% Create bivariate normal distribution of plant water delta-2-H and delta-18-O
%--------------------------------------------------------------------------
Xylem_meas = [Sap_H2_meas Sap_O18_meas]; % Measured xylem water isotopic compositions (delta-2-H and delta-18-O).
Sap_H2_mean = mean(Sap_H2_meas);
Sap_O18_mean = mean(Sap_O18_meas);
Xylem_mean = [Sap_H2_mean Sap_O18_mean]; % Mean values (delta-2-H and delta-18-O) of measured xylem water isotopic compositions.
Xylem_cov = cov(Sap_H2_meas, Sap_O18_meas); % Covariance of measured xylem water delta-2-H and delta-18-O.

% Randomly sample a normal distribution of xylem water isotopic compositions that was created using measured values.
Xylem_dist = mvnrnd(Xylem_mean, Xylem_cov, Xylem_count+1); 

%--------------------------------------------------------------------------
%                    Additional Model Settings
%--------------------------------------------------------------------------

% Set interval for computations to [0,1]
s = 3;
a = 0;
b = 1;
scale = 1;
k = 5;
n = k+s; % n is the number of intervals (spline nodes to estimate).
y = n;
omega = 1e-4; % weighting factor.

% Compute nodes (currently equi-distant)
x=zeros((n+1),1);
h = (b-a)/n;
for i=1:n+1
    x(i) = a + h*(i-1);
end

% Interpolate root distribution at specified depths
root_p=fnval(roots,x);
totalRoot=integral(@(x) fnval(roots,x),0,1);
root_pdf=fnval(roots,x)/totalRoot;

% Define depth at which soil 2H and 18O, and root water uptake pdf will be evaluated
zz=0:0.01:1; zz=zz*250; zz(1)=zz(1)+0.001; zz(end)=zz(end)-0.001; 

% Define matrices
Sap_H2_samples = zeros (Xylem_count+1,1);
Sap_O18_samples = zeros (Xylem_count+1,1);

PredQtt=zeros(Xylem_count+1,1);
Qtt_Pave1=zeros(Xylem_count+1,1);

rootwaterPDF=zeros(length(zz),1);
rootwaterPDFsum=zeros(1,Xylem_count+1);

rootwaterPDF_dist=zeros(length(zz),Xylem_count+1);

Uptake=zeros(length(zz),(Xylem_count+1));

C=zeros(Xylem_count+1,1);
D=zeros(Xylem_count+1,1);

A=zeros((s+1+6)*y,(s+1+6)*y);
rhs=zeros((s+1+6)*y,1);

SSEKernal=zeros(3,y+1);


%% ------------------------------------------------------------------------
%                           Main Program
% %------------------------------------------------------------------------

for j = 1:Xylem_count+1
    
    if j==1
        Sap_H2_samples(j) = Sap_H2_mean;
        Sap_O18_samples(j) = Sap_O18_mean;
    else
        Sap_H2_samples(j) = Xylem_dist(j,1);
        Sap_O18_samples(j) = Xylem_dist(j,2);
    end

% Add equations for the concentration and fluxes Accoding to Cook et al. (2006)

% Integrating Equation 3 with respect to time...
LeafP=-integral(@(tt) fnval(leafPsi,tt),leafStime,leafEndtime)+5.5*9.81/1000*(leafEndtime-leafStime);

% The function within the integral involving depth...
f11=@(z) cond*(LeafP+(fnval(soilPsi,z)-(Zr/100*z*9.81/1000)));

% Probability density functions
KernalType = {'NormK','BetaK','Bernstein'};
Kernal = KernalType{3};
WEIGHTS =10; % Weighting (ratio) of plant water isotopes (2-H,18-O) relative to sap flux (Qt) in the model objective function.

switch Kernal
    case 'NormK'
        band = 0.08;
        SSEKernal(1,:) = integral(@(xx) (f11(xx).*(fnval(C18O_M,xx)-...                         
            ((Sap_O18_samples+1000)*2.0052))).*normpdf(xx,x,band),0,1,'ArrayValued',true);
        SSEKernal(2,:) = integral(@(xx) (f11(xx).*(fnval(C2H_M,xx)...                           
            -((Sap_H2_samples+1000)/6.420135))).*normpdf(xx,x,band),0,1,'ArrayValued',true);
        SSEKernal(3,:)=  integral(@(xx) f11(xx).*normpdf(xx,x,band),0,1,'ArrayValued',true);    
        Aeq(1,:)=integral(@(xx) normpdf(xx,x,band),0,1,'ArrayValued',true);

    case 'BetaK'
        x(1)= x(1)+0.001; % to avoid problems with infinite values (i.e., division by zero).
        x(n+1)= x(n+1)-0.001; % to avoid problems with infinite values (i.e., division by zero).
        band = 0.23;
        kappaFun = @(x1,band,kappa) band^2*kappa.^3+(x1^3-x1^2+7*band^2*x1-...
            3*band^2)*kappa.^2+(2*x1-1)*(-x1^2+8*band^2*x1-3*band^2).*kappa...
            +band^2*(3*x1-1)*(2*x1-1)^2;
        for i=1:n+1
            kappa(i)= fsolve(@(kappa) kappaFun(x(i),band,kappa),2);
            lamda(i)= (1-x(i))/x(i)*kappa(i)+(2*x(i)-1)/x(i);
        end
        SSEKernal(1,:)= integral(@(xx) (f11(xx).*(fnval(C18O_M,xx)-...
            ((Sap_O18_samples+1000)*2.0052))).*betapdf(xx,kappa,lamda),0,1,'ArrayValued',true);
        SSEKernal(2,:)= integral(@(xx) (f11(xx).*(fnval(C2H_M,xx) -...
            ((Sap_H2_samples+1000)/6.420135))).*betapdf(xx,kappa,lamda),0,1,'ArrayValued',true);
        SSEKernal(3,:)=  integral(@(xx) f11(xx).*betapdf(xx,kappa,lamda),0,1,'ArrayValued',true');
        
        Aeq(1,:)= integral(@(xx) betapdf(xx,kappa,lamda),0,1,'ArrayValued',true');

    case 'Bernstein'
        SSEKernal(1,:) = WEIGHTS*integral(@(xx) (f11(xx).*(fnval(C18O_M,xx)-...
            ((Sap_O18_samples(j)+1000)*2.0052))).*bernstein(y, xx),0,1,'ArrayValued',true);
        SSEKernal(2,:) = WEIGHTS*integral(@(xx) (f11(xx).*(fnval(C2H_M,xx) -...
            ((Sap_H2_samples(j)+1000)/6.420135))).*bernstein(y, xx),0,1,'ArrayValued',true);
        SSEKernal(3,:) =  integral(@(xx) f11(xx).*bernstein(y, xx),0,1,'ArrayValued',true);
        for kk=1:length(x)
            SSEKernal(3+kk,:) =  bernstein(y, x(kk));
        end
        Aeq(1,:)=integral(@(xx) bernstein(y, xx),0,1,'ArrayValued',true);
        Aeq(2,:)=bernstein(y, 1); % root water uptake pdf at Zr (at bottom of the root zone).
%         Aeq(3,:)=bernstein(y, 0); % root water uptake pdf at the first point (at the surface).
end

% Equality constraints 
beq(1,1)=1;  % area under the curve = 1.
beq(2,1)=0;  % root water uptake pdf at first point (at the soil surface)= 0.
% beq(3,1)=0;  % root water uptake pdf at the last point (at the bottom of the root zone)= 0.

% Inequality for securing positive weights for Kernal density functions and coefficients for Bernstein plynomials
Alphaineq=zeros(y,y);alphabineq=zeros(y+1,1);
for i=1:y+1                     
    Alphaineq(i,i)    = -1;    
    alphabineq(i,1)   =  0;
end

%--------------------------------------------------------------------------
% Sum(Meaured-Estimated)^2 for 18-O, 2-H and Sap Flux
% Right hand side of the 18-O ,2-H, and average sap flux estimates

rhsH(1,1)=0;   % for 18-O
rhsH(2,1)=0;   % for 2-H
rhsH(3,1)=integral(@(tt)fnval(sapFlux,tt),leafStime,leafEndtime);  % For total sap flux within the day.
for kk=1:length(x)
    rhsH(3+kk,1) = root_pdf(kk);
end
 
% Convert SSE to the form: alpha'*H*alpha+f*alpha+0.5*rhsH(3,1)^2; where alpha are the unknowns.
H=SSEKernal'*SSEKernal; % SSE Matrix (of isotope and sap flux residuals).
% i.e., H =(h+H')/2;

f=-rhsH(3:(y+3),1)'*SSEKernal(3:(y+3),:); 
% i.e., f = -rhsH.*SSEKernal; % Qt(meas)*Qt(pred)
% i.e., f = f*f';

% Specify Matlab Solver
Algorithm ={'QuadProg','ConeProg'};
Prog=Algorithm{1};

switch Prog

    case 'QuadProg' % Call matlab solver quadprog.

        [sol1, fval1]=quadprog(H,f,Alphaineq,alphabineq,Aeq,beq);

        % Add the constant to the SSE
        fval=fval1+0.5*rhsH(3,1)^2;
        
    case "ConeProg"  % Call matlab solver coneprog.   % Objective function is t+f*y+0.5*C^2=fsc*[y,t]+0.5*C^2 subject.

        Asc=real(sqrtm(H(1:n+1,1:n+1)));
        Asc(n+2,n+2)=1;
        dsc=[zeros(size(f(:)));1];
        bsc=zeros(size(dsc));
        Gamma=-1;
        socConstraints=secondordercone(Asc,bsc,dsc,Gamma);

        Aineq1=[Alphaineq,zeros(n+1,1)];
        Aeq1=[Aeq, zeros(1,1)];bineq1=[alphabineq,zeros(n+1,1)];
        fsc=[f,1]; % t=1
        options=optimoptions('Coneprog','Display','iter','LinearSolver','normal');
        [sol2,fval]=coneprog(fsc,socConstraints,Aineq1,alphabineq,Aeq1,beq,[],[],options);
        sol1=sol2(1:end-1);        
        fval1=fval+0.5*rhsH(3:(y+3),1)'*rhsH(3:(y+3),1)+0.5;
end

%--------------------------------------------------------------------------

% Total xylem flux density
Qtt_Mave1=rhsH(3,1)/(leafEndtime-leafStime); % Mean measured xylem flux.

PredQtt(j)=SSEKernal(3,1:y+1)*sol1; 
Qtt_Pave1(j)=PredQtt(j)/(leafEndtime-leafStime); % Predicted xylem flux.

%% ------------------------------------------------------------------------
%   Calculate the root water uptake pdf and store depth in 'zz' and pdf value in 'rootwaterPDF'
%  ------------------------------------------------------------------------

for kk=1:length(zz)
    sum=0;
    if strcmp(Kernal,'BetaK')==1
        sum = betapdf(zz(kk)/250,kappa,lamda)*sol1/250;
    elseif strcmp(Kernal,'NormK')==1
        sum= arrayfun(@(xx1) normpdf(zz(kk)/250,xx1,band),x)'*sol1/250;
    elseif strcmp(Kernal,'Bernstein')==1
        bernsteinK=bernstein(y,zz(kk)/250);
        sum =bernsteinK*sol1/250;
    end
    
    rootwaterPDF(kk)=sum;
    
end
clear sum

% Normalize water uptake pattern, to ensure unity of PDF
rootwaterPDFsum(j)=sum(rootwaterPDF);
rootwaterPDF_dist(:,j)=rootwaterPDF;
rootwaterPDF_dist(:,j)=rootwaterPDF_dist(:,j)./rootwaterPDFsum(1,j);


% Calculate predicted xylem 18-O and 2-H 
PredO_H = zeros(2,y+1);
if strcmp(Kernal,'BetaK')==1
    PredO_H(1,:) = integral(@(xx) f11(xx).*(fnval(C18O_M,xx)).*...
        betapdf(xx,kappa,lamda),0,1,'ArrayValued',true)/PredQtt;
    PredO_H(2,:) = integral(@(xx) f11(xx).*fnval(C2H_M,xx).*...
        betapdf(xx,kappa,lamda),0,1,'ArrayValued',true)/PredQtt;
elseif strcmp(Kernal,'NormK')==1
    PredO_H(1,:) = integral(@(xx) f11(xx).*(fnval(C18O_M,xx)).*...
        normpdf(xx,x,band),0,1,'ArrayValued',true)/PredQtt;
    PredO_H(2,:) = integral(@(xx) f11(xx).*fnval(C2H_M,xx).*...
        normpdf(xx,x,band),0,1,'ArrayValued',true)/PredQtt;
elseif strcmp(Kernal,'Bernstein')==1
    PredO_H(1,:) = integral(@(xx) f11(xx).*(fnval(C18O_M,xx)).*...
        bernstein(y, xx),0,1,'ArrayValued',true)/PredQtt(j);
    PredO_H(2,:) = integral(@(xx) f11(xx).*fnval(C2H_M,xx).*...
        bernstein(y, xx),0,1,'ArrayValued',true)/PredQtt(j);        
end

% Convert predicted xylem water 2-H and 18-O compositions from units of ppm to per mil
C(j)=PredO_H(2,:)*sol1*6.420135-1000; % predicted xylem water 2-H.
D(j)=PredO_H(1,:)*sol1/2.0052-1000; % predicted xylem water 18-O.

end

%% --------------------------------------------------------------------
%                       Posterior Analysis
% %--------------------------------------------------------------------

% Mean root water uptake
rootwaterPDF_dist_mean = mean(rootwaterPDF_dist,2); 
rootwaterPDF_dist_mean_profile = [zz' rootwaterPDF_dist_mean]; 
rootwaterPDF_dist_mean_profile_curve=csape(rootwaterPDF_dist_mean_profile(:,1)/Zr,[0,rootwaterPDF_dist_mean_profile(:,2)',0],[2 2]);

% Uncertainty of root water uptake
rootwaterPDF_dist_stddev=std(rootwaterPDF_dist,0,2);
quantiles=quantile(rootwaterPDF_dist,[0.025 0.1587 0.25 0.5 0.75 0.8413 0.975],2);


Predicted=[Sap_H2_samples Sap_O18_samples C D Qtt_Pave1];
RWU_meas=[zz',rootwaterPDF_dist(:,1)]; % Predicted water uptake profile obtained using the mean of the measured 2-H and 18-O plant water isotopic compositions.
RWU_means = [RWU_meas,rootwaterPDF_dist_mean];  % Predicted water uptake profile obtained using the bivariate distribution of 2-H and 18-O generated using measured plant water isotopic compositions (i.e., MC simulation)


%% ------------------------------------------------------------------------
%                     Plotting Results
% %------------------------------------------------------------------------
fig=figure ('Name','Plots');
ColorChoice=[0.3467,    0.5360,    0.6907;...
    0.9153,    0.2816,    0.2878;...
    0.4416,    0.7490,    0.4322;...
    1.0000,    0.5984,    0.2000;...
    0.6769,    0.4447,    0.7114;...
    0     ,    0     ,    0     ];
LinSty={'-','--',':','-.'};
numLineStyle=3;
tiledlayout(1,8);
%--------------------------------------------------------------------------
ax1=nexttile; % Measured soil water potential.
plot(SoilPsi(:,2),SoilPsi(:,1),'k:d',fnval(soilPsi,psidepth/Zr),psidepth,'k-','MarkerSize',6,'MarkerFaceColor','k','LineWidth',2);
set(gca, 'YDir','reverse')
title(ax1,'Soil Potential')
ylabel(ax1,'Depth (cm)')
xlabel(ax1,'Matric Potential (MPa)')
axis(ax1,[-inf 0 0 250])
box on
grid on
legend('Measured','Interpolated')
legend('location','south')
%--------------------------------------------------------------------------
ax2=nexttile; % Measured leaf water potential.
LeafPotentialTimeS=timeseries(LeafPsi(:,2),LeafPsi(:,1),"Name",'Measured');
LeafPotentialTimeS.Name='Leaf \psi (MPa)';
LeafPotentialTimeS.TimeInfo.Units='Hours';
plot(LeafPotentialTimeS,'k:s','MarkerSize',6,'MarkerFaceColor','k','LineWidth',2);
hold on
LeafPotentialTimeS1=timeseries(fnval(leafPsi,leafPsitime), leafPsitime,"Name",'interpolated');
LeafPotentialTimeS1.Name='Leaf \psi (MPa)';
LeafPotentialTimeS1.TimeInfo.Units='Hours';
plot(LeafPotentialTimeS1,'k-','MarkerSize',6,'MarkerFaceColor','k','LineWidth',2);
title(ax2,'Leaf Potential')
ylabel(ax2,'Potential (Mpa)')
xlabel(ax2,'Time (hours)')
box on
grid on
legend('Measured','Interpolated')
legend('location','south')
%--------------------------------------------------------------------------
ax3=nexttile; % Measured root distribution.
hold on
plot(Roots(:,2),Roots(:,1),'k:d',fnval(roots,psidepth/Zr),psidepth,'k-','MarkerSize',6,'MarkerFaceColor','k','LineWidth',2);
set(gca,'YDir','reverse');
ylim([0 Zr])
title(ax3,'Root Distribution')
xlabel(ax3,'Root Length Density')
ylabel(ax3,'Depth(cm)')
box on
grid on
legend('Measured','Interpolated')
legend('location','south')
%--------------------------------------------------------------------------
ax4=nexttile; % Measured and predicted sap flux density.
k=1;
SapFluxTimeS=timeseries(Qtt_Pave1(j),leafStime:0.2:leafEndtime,"Name",'Predicted average');
SapFluxTimeS.Name='Sap Flux Density(cm h^-^1)';
SapFluxTimeS.TimeInfo.Units='Hours';
plot(SapFluxTimeS, 'color', 'r' );
hold on
Qtt_Mave1=rhsH(3,1)/(leafEndtime-leafStime);
k=k+1;
SapFluxTimeS=timeseries(Qtt_Mave1,leafStime:0.2:leafEndtime,"Name",'Measured average');
SapFluxTimeS.Name='Sap Flux Density(cm h^-^1)';
SapFluxTimeS.TimeInfo.Units='Hours';
plot(SapFluxTimeS, 'color', 'b');
SapFluxTimeS=timeseries(SapFlux(:,2),SapFlux(:,1),"Name",'Measured');
SapFluxTimeS.Name='Sap Flux Density(cm h^-^1)';
SapFluxTimeS.TimeInfo.Units='Hours';
p1=plot(SapFluxTimeS,'k:v','MarkerSize',6,'MarkerFaceColor','k','LineWidth',2);
title(ax4,'Sap Flux')
ylabel(ax4,'Sap Flux Density (cm h^-^1)')
xlabel(ax4,'Time (hours)')
grid on
legendinfo{1}="Predicted Mean";
legendinfo{2}="Measured Mean";
legendinfo{3}="Measured";
legend(legendinfo,'location','south')
box on
grid on
hold off
%--------------------------------------------------------------------------
ax5=nexttile; % Soil water 2-H composition as a function depth and the corresponding plant water delta-2-H.
hold on
plot(ax5,Delta(:,2),Delta(:,1),'k:*','LineWidth',2);
plot(ax5,Delta(1,4)*ones(length(psidepth),1),psidepth,'g:pentagram','MarkerSize',4,'MarkerIndices',1:30:length(psidepth),'LineWidth',2);
plot(C(1,1).*ones((length(psidepth)),1),psidepth,'color','r');
set(gca,'YDir','reverse')
legendinfo{1}='Soil Water';
legendinfo{2}='Measured Xylem';
legendinfo{3}='Predicted Xylem';
legend(legendinfo,'location','south')
grid on
title(ax5,'\delta^2H')
ylabel(ax5,'Depth (cm)')
xlabel(ax5,['\delta^2H(',char(8240),')'])
box on
grid on
hold off
%--------------------------------------------------------------------------
ax6=nexttile; % Soil water 18-O composition as a function depth and the corresponding plant water delta-18-O.
hold on
plot(ax6,Delta(:,3),Delta(:,1),'k:*','LineWidth', 2);
plot(ax6,Delta(1,5)*ones(length(psidepth),1),psidepth,'g:pentagram','MarkerSize',4,'MarkerIndices',1:30:length(psidepth),'LineWidth',2);
plot(D(1,1).*ones(length(psidepth),1),psidepth, 'color', 'r'); 
set(gca,'YDir','reverse')
legendinfo{1}='Soil Water';
legendinfo{2}='Measured Xylem';
legendinfo{3}='Predicted Xylem';
legend(legendinfo,'location','south')
title(ax6,'\delta^1^8O')
ylabel(ax6,'Depth (cm)')
xlabel(ax6,['\delta^1^8O (', char(8240),')'])
box on
grid on
hold off
%--------------------------------------------------------------------------
ax7=nexttile; % Measured root distribution, predicted mean water uptake, and predicted 2.5% and 97.5% quantiles of water uptake.
hold on
root_pdfz=fnval(roots,zz/Zr)/totalRoot; % Measured root distribution.
plot(ax7,root_pdfz/Zr,zz,'k-','lineWidth',2);
plot(ax7,RWU_meas(:,2),zz,'r-','LineWidth',2);
plot(ax7,rootwaterPDF_dist_mean,zz,'r--','LineWidth',2);
plot(ax7,quantiles(:,1),zz,'b:','lineWidth',2); % 2.5% quantile of predicted water uptake.
plot(ax7,quantiles(:,7),zz,'b:','lineWidth',2); % 97.5% quantile of predicted water uptake.
set(gca,'YDir','reverse');
ylim([0 Zr])
legend({'Measured Roots','Measured RWU','Predicted RWU','2.5%Q','97.5%Q'},'location','south')
title(ax7,'Water Uptake Pattern')
xlabel(ax7,'Probability Density (cm^-^1)')
ylabel(ax7,'Depth(cm)')
box on
grid on
hold off
%--------------------------------------------------------------------------
ax8=nexttile; % Mean predicted RWU (from MC simulation)
hold on
plot(rootwaterPDF_dist_mean_profile(:,2),rootwaterPDF_dist_mean_profile(:,1),'k:d',fnval(rootwaterPDF_dist_mean_profile_curve,psidepth/Zr),psidepth,'k-','MarkerSize',6,'MarkerFaceColor','k','LineWidth',2);
set(gca,'YDir','reverse');
ylim([0 Zr])
title(ax3,'Mean Predicted RWU')
xlabel(ax3,'Water Uptake (%)')
ylabel(ax3,'Depth(cm)')
box on
grid on
legend('Measured','Interpolated')
legend('location','south')

%% ------------------------------------------------------------------------
%                 Write output data to excel files
% %------------------------------------------------------------------------

header = {'Depth (cm)','Meas','Pred'};
RWU_means = [header; num2cell(RWU_means)];

header_quantiles = {'2.5%','-Std Dev','1st quartile','Median','3rd quartile','+Std Dev','97.5%'};
quantiles = [header_quantiles; num2cell(quantiles)];

header = {'H2_samples','O18_samples','H2_pred','O18_pred','Qt_pred'};
Predicted = [header; num2cell(Predicted)];

writecell(RWU_means, [fileName(1:end-5),'_RWU_means.xlsx']);
writecell(quantiles, [fileName(1:end-5),'_RWU_quantiles.xlsx']);
writecell(Predicted, [fileName(1:end-5),'_predicted_values.xlsx']);

%% ------------------------------------------------------------------------
%                       Model Functions
% %------------------------------------------------------------------------

function Fr= Beta(z,kappa,lamda)
if(z<=0 ||z>=1)
    Fr=0;
else
    lnFr = (kappa-1)*log(z)+(lamda-1)*log(1-z)-betaln(kappa,lamda);
    Fr = exp(lnFr);
end
end

function [c ceq]=nonlcon(alpha)
n=length(alpha)/10;
for i=1:n
    c(i)=-alpha(7*n+i)*alpha(9*n+i)+alpha(8*n+i)^2;   % -xi*zi+yi^2<=0
    c(n+i)=-alpha(4*n+i)*alpha(6*n+i)+alpha(5*n+i)^2; % -si*wi+vi^2<=0
end
ceq=[];
end

function bern = bernstein(n,x)
%
%% BERNSTEIN evaluates the Bernstein polynomials defined on [0,1].
% 
% The Bernstein polynomials are assumed to be based on [0,1].
% 
%  Formula:
%    B(N,I)(X) = [N!/(I!*(N-I)!)] * (1-X)^(N-I) * X^I
%
%%      First values:
%% 
%       B(0,0)(X) = 1
% 
%%
%       B(1,0)(X) =      1-X
%       B(1,1)(X) =               X
%
%%
%       B(2,0)(X) =     (1-X)^2
%       B(2,1)(X) = 2 * (1-X)   * X
%       B(2,2)(X) =               X^2
%
%%
%       B(3,0)(X) =     (1-X)^3
%       B(3,1)(X) = 3 * (1-X)^2 * X
%       B(3,2)(X) = 3 * (1-X)   * X^2
%       B(3,3)(X) =               X^3
%
%%
% 
%       B(4,0)(X) =     (1-X)^4
%       B(4,1)(X) = 4 * (1-X)^3 * X
%       B(4,2)(X) = 6 * (1-X)^2 * X^2
%       B(4,3)(X) = 4 * (1-X)   * X^3
%       B(4,4)(X) =               X^4
%
%%   Special values:
%%
%    B(N,I)(X) has a unique maximum value at X = I/N.
%
%%
%    B(N,I)(X) has an I-fold zero at 0 and and N-I fold zero at 1.
%
%%
%    B(N,I)(1/2) = C(N,K) / 2^N
%
%%
%    For a fixed X and N, the polynomials add up to 1:
%    Sum ( 0 <= I <= N ) B(N,I)(X) = 1
%
%%   Licensing:
%%   
%    This code is distributed under the GNU LGPL license.
%
%%
%    Modified:
%    22 July 2004
%    Author: John Burkardt
%
%%   Parameters:
%%
%    Input, integer N, the degree of the Bernstein polynomials to be used.  
%    For any N, there is a set of N+1 Bernstein polynomials,
%    each of degree N, which form a basis for polynomials on [0,1].
%
%    Input, real X, the evaluation point.
%
%    Output, real BERN(1:N+1), the values of the N+1 Bernstein polynomials at X.
%
%%

if ( n == 0 )
    
    bern(1) = 1.0;
    
elseif ( 0 < n )
    
    bern(1) = 1.0 - x;
    bern(2) = x;
    
    for i = 2 : n
        bern(i+1) = x * bern(i);
        for j = i-1 : -1 : 1
            bern(j+1) = x * bern(j) + ( 1.0 - x ) * bern(j+1);
        end
        bern(1) = ( 1.0 - x ) * bern(1);
    end
    
end
end