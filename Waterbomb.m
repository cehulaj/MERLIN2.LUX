%% ========================== Waterbomb ================================ %%
clear all; close all;  clc;

%% Define geomtry
a = 100;
b = 70;
h = 20;

Node = [ 0,  0,  h;
         b,  0,  0;
         0,  b,  0;
        -b,  0,  0;
         0, -b,  0;
         a,  a,  0;
        -a,  a,  0;
        -a, -a,  0;
         a, -a,  0;];
  
Panel = {[1,2,6];[1,6,3];[1,3,7];[1,7,4];[1,4,8];[1,8,5];[1,5,9];[1,9,2]};

%% Interesting material parameters
% Force mode
%FDef = 15; EY = 3e3; t = 0.27; Ls = 50; EY2 = 3e3; Ls2 = Ls; D_lmd = 0.001; interv = 20; %Paper
FDef = 0.2; EY = 10; t = 1; Ls = 50; EY2 = 30; Ls2 = Ls; D_lmd = 0.01; interv = 20; %LASMP

%% Set up boundary conditions
m = size(Node,1);
Supp = [ 1, 1, 1, 0;
         6, 0, 0, 1;
         7, 0, 0, 1;
         8, 0, 0, 1;
         9, 0, 0, 1];

% Force mode
Load = [ 1, 0, 0, -FDef];
% Displacement mode
%Load = [ 1, 0, 0, -150]; interv = 4; % LASMP

%% Adopt generalized N5B8 model
% Force mode

AnalyInputOpt = struct(...
    'ModelType','N5B8',...
    'MaterCalib','auto',...
    'ModElastic', EY,...
    'Poisson', 0.33,...
    'Thickness', t,... 
    'LScaleFactor', Ls,...
    'ModElastic2', EY2,...
    'LoadType','Force',...
    'InitialLoadFactor', D_lmd,...
    'MaxIcr', 5000,...
    'StopCriterion',@(Node,U,icrm,lmd)(abs(lmd)>2));
%}
% Displacement mode
%{
AnalyInputOpt = struct(...
    'ModelType','N5B8',...
    'MaterCalib','auto',...  
    'ModElastic', 10,...
    'Poisson', 0.33,...
    'Thickness', 1,... 
    'LScaleFactor', 50,...
    'ModElastic2', 30,...
    'LScaleFactor2', 50,...
    'LoadType','Displacement',...  % Displacement load
    'InitialLoadFactor', 0.004,...
    'MaxIcr', 1000,...
    'DispStep',260);
%}

%% Perform analysis
% Assemble input data
[truss, angles, AnalyInputOpt] = PrepareData(Node,Panel,Supp,Load,AnalyInputOpt);

% Specify initial deformation state
truss.U0 = zeros(3*size(truss.Node,1),1);

% Perform path-following analysis using MGDCM
[Uhis,Fhis,angles,truss,IraIcrm] = PathAnalysis(truss,angles,AnalyInputOpt);
% After iraicrm, the origami is iradiated, before, not iradiated

% Clean output data
Uhis = real(Uhis);
Fhis = real(Fhis);
STAT = PostProcess(Uhis,truss,angles,IraIcrm); 

%% Visualize simulation
instdof = [1, 3]; endicrm = size(Uhis,2);

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

%  Force vs RealH
    subplot(2,2,4);
    RelH = sign(instdof(2))*Uhis((instdof(1)*3-(3-abs(instdof(2)))),:)/h+1;
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
    title('Waterbomb','fontsize',14,'fontweight','normal');
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
    