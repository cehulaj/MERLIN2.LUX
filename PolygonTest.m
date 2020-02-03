
%% ========================== Polygon test ============================== %%
clear all; close all; clc;

%% Define geometry
Node = [ 0,  0,  0;
         0, 10,  0;
        12, 15,  0;
        30, 10,  0;
        30,  0,  0;
        18, -5,  0]*2;

Panel = {[1:6]};

%% Set up boundary conditions
m = size(Node,1);
Supp = [ 1, 1, 1, 1;
         2, 1, 0, 1;
         4, 0, 0, 1];
Load = [ 5, 0, 0,  -0.5];

%% Adopt generalized N5B8 model
% Auto mode

AnalyInputOpt = struct(...
    'ModelType','N5B8',...
    'MaterCalib','auto',...
    'ModElastic', 10,...
    'Poisson', 0.33,...
    'Thickness', 1,... 
    'LScaleFactor', 50,...
    'ModElastic2', 30,...
    'LScaleFactor2', 50,...
    'LoadType','Force',...
    'InitialLoadFactor', 0.0000005,...
    'MaxIcr', 2000,...
    'StopCriterion',@(Node,U,icrm,lmd)(abs(lmd)>2));
%}
%% Adopt generalized N4B5 model
% Manual mode
%{
AnalyInputOpt = struct(...
    'ModelType','N4B5',...
    'MaterCalib','manual',...
    'BarCM', @(Ex)Ogden(Ex, 5e3),...
    'BarCM2', @(Ex2)Ogden(Ex2, 5e3),...
    'Abar', 2.4652,...
    'Kb',0.9726*([38.4187;38.4187;41.7612]/0.127).^(1/3)./[38.4187;38.4187;41.7612],...
    'Kf',1,...
    'Kb2',0.9726*([38.4187;38.4187;41.7612]/0.127).^(1/3)./[38.4187;38.4187;41.7612],...
    'Kf2',1,...
    'RotSprBend', @SuperLinearBend,...
    'RotSprFold', @(he,h0,Kf,L0)EnhancedLinear(he,h0,Kf,L0,15,345),...
    'RotSprBend2', @SuperLinearBend,...
    'RotSprFold2', @(he,h02,Kf2,L0)EnhancedLinear(he,h02,Kf2,L0,15,345),...
    'LoadType','Force',...
    'InitialLoadFactor', 0.000005,...
    'MaxIcr', 1000,...
    'StopCriterion',@(Node,U,icrm,lmd)(abs(lmd)>10));
%}

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
instdof = [5,-3];
interv = 5; endicrm = size(Uhis,2);

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

%  Force vs Displacement
    subplot(2,2,4);
    dsp = sign(instdof(2))*Uhis((instdof(1)*3-(3-abs(instdof(2)))),:);
    plot(dsp(1:IraIcrm),Fhis(1:IraIcrm),'-','Color','k','LineWidth',1.5)
    hold on
    plot(dsp(IraIcrm:end),Fhis(IraIcrm:end),'--','Color','k','LineWidth',1.5)
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
    %rotate3d on
    title('Polygon','fontsize',14,'fontweight','normal');
    xlabel('x (mm)','fontsize',12);
    ylabel('y (mm)','fontsize',12);
    zlabel('z (mm)','fontsize',12);   

%% Plot diagrams
%  Force vs Displacement
    figure()
    plot(dsp(1:IraIcrm),Fhis(1:IraIcrm),'Color',[0 0 1],'LineWidth',1.5)
    hold on
    plot(dsp(IraIcrm:end),Fhis(IraIcrm:end),'Color',[0.3010 0.7450 0.9330],'LineWidth',1.5)
    axis tight
    xlabel('Displacement','fontsize',14)
    ylabel('Load Factor','fontsize',14);
    
%  Strain energy vs Displacement
    figure()
    plot(dsp(1:IraIcrm),STAT.PE(1:IraIcrm),'Color',[0.6350 0.0780 0.1840],'LineWidth',1.5)
    hold on
    plot(dsp(IraIcrm:end),STAT.PE(IraIcrm:end),'Color',[0.9290 0.6940 0.1250],'LineWidth',1.5)
    axis tight
    xlabel('Displacement','fontsize',14)
    ylabel('Stored Energy','fontsize',14);
