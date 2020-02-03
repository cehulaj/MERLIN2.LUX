%% ============================ BeelzeN ================================ %%
clear all;  close all;  clc;

%% Interesting material parameters
%FDef = 300; EY = 10; EY2 = 30; t = 1; Ls = 50; Ls2 = 50; D_lmd = 0.0000007; %1
%FDef = 250; EY = 10; EY2 = 30; t = 1; Ls = 50; Ls2 = 50; D_lmd = 0.0000025; %2
%FDef = 200; EY = 10; EY2 = 30; t = 1; Ls = 50; Ls2 = 50; D_lmd = 0.0000003; %3
FDef = 2; EY = 10; EY2 = 30; t = 1; Ls = 50; Ls2 = 50; D_lmd = 0.00001; %3

%% Define geometry
N = 23;
%R1 = 70; R2 = 60; R3 = 80; %1
%R1 = 70; R2 = 70; R3 = 95; %2
R1 = 70; R2 = 60; R3 = 75; %3

[Node, Panel] = GetBeelzeN(N,R1,R2,R3);

%% Set up boundary conditions
Supp = ones(N+1,4);
for n = 1:N
    Supp(n,1) = N+n;
end
Supp(N+1,1) = 3*N+1;

Load = [ 1, 0, 0, 0.1];

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
    'AdaptiveLoad',@LoadFun,...
    'OrigamiType', 'BeelzeN',...
    'N',N,...
    'FDef', FDef,...
    'InitialLoadFactor', D_lmd,...
    'MaxIcr', 2000,...
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
% First node: z component
%instdof = [1, 3];
interv = 10; endicrm = size(Uhis,2);

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
   % view(-2,16)
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

%  Force vs DTh = Th - Th01
    subplot(2,2,4);
    DTh = 180/pi*acos((Uhis(1,:)+R2)./sqrt((R2+Uhis(1,:)).^2+Uhis(3,:).^2));
    plot(DTh(1:IraIcrm),Fhis(1:IraIcrm),'-','Color','k','LineWidth',1.5)
    hold on
    plot(DTh(IraIcrm:end),Fhis(IraIcrm:end),'--','Color','k','LineWidth',1.5)
    axis tight
    xlabel('\theta - \theta_{01} (°)','fontsize',12)
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
    title('BeelzeN','fontsize',14,'fontweight','normal');
    xlabel('x (mm)','fontsize',12);
    ylabel('y (mm)','fontsize',12);
    zlabel('z (mm)','fontsize',12);
    
%% Plot diagrams
%  Force vs DTh = Th - Th01
    figure()
    plot(DTh(1:IraIcrm),Fhis(1:IraIcrm),'Color',[0 0 1],'LineWidth',1.5)
    hold on
    plot(DTh(IraIcrm:end),Fhis(IraIcrm:end),'Color',[0.3010 0.7450 0.9330],'LineWidth',1.5)
    axis tight
    xlabel('Rotational Angle','fontsize',14)
    ylabel('Load Factor','fontsize',14);
    
%  Strain energy vs DTh = Th - Th01
    figure()
    plot(DTh(1:IraIcrm),STAT.PE(1:IraIcrm),'Color',[0.6350 0.0780 0.1840],'LineWidth',1.5)
    hold on
    plot(DTh(IraIcrm:end),STAT.PE(IraIcrm:end),'Color',[0.9290 0.6940 0.1250],'LineWidth',1.5)
    axis tight
    xlabel('Rotational Angle','fontsize',14)
    ylabel('Stored Energy','fontsize',14);
