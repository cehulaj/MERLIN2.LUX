%% ============================= Kresling ============================== %%
clear all;  close all;  clc;

%% Define geometry
N = 12; R = 50; h = 50; lyr = 4; phi = 2*pi/N; H0 = (lyr-1)*h; %1
%N = 12; R = 50; h = 50; lyr = 6; phi = 2*pi/N; H0 = (lyr-1)*h; %2
%N =  8; R = 50; h = 50; lyr = 6; phi = 2*pi/N; H0 = (lyr-1)*h; %3
% N-gon prism; R - radius of the pattern, h - height of each layer; 
% lyr - number of inter-layer planes; phi - twisting angle of the prism in each layer;
% H0 - height of the pattern

[Node, Panel] = GetDiSym(N,R,h,lyr,phi);

%% Interesting material parameters
FDef = 1.5; EY = 10; t = 1.0; Ls = 5; EY2 = 50; Ls2 = Ls; D_lmd = 0.002; MaxIcr = 20000; interv = 50; %1
%FDef = 1.5; EY = 10; t = 1.0; Ls = 50; EY2 = 30; Ls2 = Ls; D_lmd = 0.002; MaxIcr = 10000; interv = 50; %1
%FDef = 1.5; EY = 10; t = 1.0; Ls = 5; EY2 = 30; Ls2 = Ls; D_lmd = 0.02; MaxIcr = 5000; interv = 5; %1
%FDef = 1.5; EY = 10; t = 1.0; Ls = 5; EY2 = 30; Ls2 = Ls; D_lmd = 0.01; MaxIcr = 10000; interv = 10; %2
%FDef = 0.5; EY = 10; t = 1.0; Ls = 5; EY2 = 30; Ls2 = Ls; D_lmd = 0.01; MaxIcr = 10000; interv = 10; %3

%% Set boundary condition
indsupp = find(Node(:,3)<0.01);
nsupp = numel(indsupp);
Supp = [ indsupp(1), 1, 1, 1;
         indsupp(2), 1, 1, 1;
         indsupp(3:end), zeros(nsupp-2,1)+1, zeros(nsupp-2,1)+1, ones(nsupp-2,1);];
m = size(Node,1);
indp = find(abs(Node(:,3)-max(Node(:,3)))<1e-5); %indp = indp(3);
npp = numel(indp);
Load = [indp, 0*ones(npp,1), 0*ones(npp,1), -FDef*ones(npp,1);];

%% Adopt generalized N5B8 model
AnalyInputOpt = struct(...
    'ModelType','N5B8',...
    'MaterCalib','auto',...
    'ModElastic', EY,...
    'Poisson', 0.33,...
    'Thickness', t,... 
    'LScaleFactor', Ls,...
    'ModElastic2', EY2,...
    'LScaleFactor2', Ls2,...
    'LoadType','Force',...
    'InitialLoadFactor', D_lmd,...
    'MaxIcr', MaxIcr,...
    'StopCriterion',@(Node,U,icrm,lmd)(abs(lmd)>2));
 
%% Perform analysis
% Assemble input data
[truss, angles, AnalyInputOpt] = PrepareData(Node,Panel,Supp,Load,AnalyInputOpt);

% Specify initial deformation state
truss.U0 = zeros(3*size(truss.Node,1),1);
     
% Perform path-following analysis using MGDCM
[Uhis,Fhis,angles,truss,IraIcrm] = PathAnalysis(truss,angles,AnalyInputOpt);
% After IraIcrm, the origami is iradiated, before, not iradiated

% Clean output data
Uhis = real(Uhis);
Fhis = real(Fhis);
STAT = PostProcess(Uhis,truss,angles,IraIcrm); 

%% Visualize simulation
instdof = [lyr*N, 3]; endicrm = size(Uhis,2);

%VisualFold(Uhis(:,1:interv:endicrm),truss,angles,Fhis(1:interv:endicrm,:),instdof)

% If do not need load-displacement diagram:
VisualFold(Uhis(:,1:interv:endicrm),truss,angles,[],[])

% To visualize the strains in bars:
%VisualFold(Uhis(:,1:10:endicrm),truss,angles,[],[],'IntensityMap','Edge','IntensityData',STAT.bar.Sx(:,1:10:endicrm),'ShowInitial','off')

%% Plot configurations and load displacement diagram
% Plot initial configutaion
   figure()
   
   subplot(2,2,1);
   PlotOri(truss.Node,angles.Panel,truss.Trigl,'PanelColor',[0.9 0.9 0.9]);
   axis equal; axis off;
   axis tight
   title('Original shape','fontsize',12,'fontweight','normal')
   camproj('perspective')
   light
   view(-2,16)
   %rotate3d on
   
% Plot irradiated configuration
   Ui = Uhis(:,IraIcrm);
   Nodei = truss.Node; 
   Nodei(:,1) = truss.Node(:,1)+Ui(1:3:end);
   Nodei(:,2) = truss.Node(:,2)+Ui(2:3:end);
   Nodei(:,3) = truss.Node(:,3)+Ui(3:3:end);
   subplot(2,2,2);
   PlotOri(Nodei,angles.Panel,truss.Trigl,'PanelColor',[0.9 0.9 0.9]);
   axis equal; axis off;
   axis tight
   title('Irradiated shape','fontsize',12,'fontweight','normal')
   camproj('perspective')
   light
   view(-2,16)
   %rotate3d on
   
% Plot final configuration
   Ux = Uhis(:,end);
   Nodew = truss.Node;
   Nodew(:,1) = truss.Node(:,1)+Ux(1:3:end);
   Nodew(:,2) = truss.Node(:,2)+Ux(2:3:end);
   Nodew(:,3) = truss.Node(:,3)+Ux(3:3:end); 
   subplot(2,2,3);
   PlotOri(Nodew,angles.Panel,truss.Trigl,'PanelColor',[0.9 0.9 0.9]);
   axis equal; axis off;
   axis tight
   title('Temporary shape','fontsize',12,'fontweight','normal')
   camproj('perspective')
   light
   view(-2,16)
   %rotate3d on

%  Force vs RelH
    subplot(2,2,4);
    RelH = sign(instdof(2))*Uhis((instdof(1)*3-(3-abs(instdof(2)))),:)/H0+1;
    plot(RelH(1:IraIcrm),Fhis(1:IraIcrm),'-','Color','k','LineWidth',1.5)
    hold on
    plot(RelH(IraIcrm:end),Fhis(IraIcrm:end),'--','Color','k','LineWidth',1.5)
    axis tight
    xlabel('\eta','fontsize',12)
    ylabel('\epsilon','fontsize',12);

%% Plot representational scheme
    figure()
    PlotOri(truss.Node,angles.Panel,truss.Trigl,'ShowNumber','on');
    axis equal; 
    %axis off;
    axis tight
    camproj('perspective')
    %light
    view(-2,16)
    rotate3d on
    title('Kresling','fontsize',14,'fontweight','normal');
    xlabel('x (mm)','fontsize',12);
    ylabel('y (mm)','fontsize',12);
    zlabel('z (mm)','fontsize',12);
    
%% Plot diagrams   
%  Force vs RelH
    figure()
    plot(RelH(1:IraIcrm),Fhis(1:IraIcrm),'Color',[0 0 1],'LineWidth',1.5)
    hold on
    plot(RelH(IraIcrm:end),Fhis(IraIcrm:end),'Color',[0.3010 0.7450 0.9330],'LineWidth',1.5)
    axis tight
    xlabel('Relative Height','fontsize',14)
    ylabel('Load Factor','fontsize',14);
    
%  Strain energy vs RelH
    figure()
    plot(RelH(1:IraIcrm),STAT.PE(1:IraIcrm),'Color',[0.6350 0.0780 0.1840],'LineWidth',1.5)
    hold on
    plot(RelH(IraIcrm:end),STAT.PE(IraIcrm:end),'Color',[0.9290 0.6940 0.1250],'LineWidth',1.5)
    axis tight
    xlabel('Relative Height','fontsize',14)
    ylabel('Stored Energy','fontsize',14);
