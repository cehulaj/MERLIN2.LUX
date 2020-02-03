%% ======================= Miura bending =============================== %%
clear all; close all; clc;

%% Define geometry
% Import geometry in OBJ format
[Node, Panel] = ReadOBJ('GMiura_FreeformOri.obj');

%% Set up boundary conditions
m = size(Node,1);
Supp = [     1, 0, 1, 1;
             4, 1, 1, 1;
             8, 0, 1, 1;
            22, 0, 0, 1;
            59, 0, 0, 0;
            63, 1, 0, 1;
            81, 1, 0, 1;];
indp = 59;   
% Force mode
%ff = -3*ones(length(indp),1);
%Load = [indp,-3*ff/5,ff,zeros(length(indp),1)];
% Displacement mode
ff = -30*ones(length(indp),1);
Load = [indp,-1*ff/5,ff,zeros(length(indp),1)];

%% Adopt generalized N5B8 model
% Displacement mode

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
    'InitialLoadFactor', 0.002,...
    'MaxIcr', 1500,...
    'DispStep',200);
%}
% Force mode
%{
AnalyInputOpt = struct(...
    'ModelType','N5B8',...
    'MaterCalib','auto',...  
    'ModElastic', 1e3,...
    'Poisson', 0.33,...
    'Thickness', 0.37,... 
    'LScaleFactor', 50,...
    'ModElastic2', 1e3,...
    'LoadType','Force',...
    'InitialLoadFactor', 0.0001,...
    'MaxIcr', 1500,...
    'StopCriterion',@(Node,U,icrm,lmd)(abs(lmd)>2));
%}

%% Adopt generalized N4B5 model
%{
AnalyInputOpt = struct(...
    'ModelType','N4B5',...
    'MaterCalib','manual',... 
    'BarCM', @(Ex)Ogden(Ex, 5e3),...
    'Abar', 0.1,...
    'Kb',1,...
    'Kf',0.1,...
    'RotSprBend', @(he,h0,Kb,L0)EnhancedLinear(he,h0,Kb,L0,15,345),...
    'RotSprFold', @(he,h0,Kf,L0)EnhancedLinear(he,h0,Kf,L0,15,345),...
    'BarCM2', @(Ex2)Ogden(Ex2, 5e3),...
    'Kb2',1,...
    'Kf2',0.1,...
    'LoadType','Force',...
    'InitialLoadFactor', 0.005,...
    'MaxIcr', 1500,...
    'StopCriterion',@(Node,U,icrm,lmd)(abs(lmd)>2));
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
instdof = [indp(1),-2];
interv = 10; endicrm = size(Uhis,2);
% Animation monitoring node-wise change
VIntensityDataInten = zeros(size(truss.Node,1),size(Uhis,2));
IntensityDataM = bsxfun(@times,STAT.bar.Sx,truss.A);
for k = 1:size(Uhis,2)
    IntensityDataIntenk = sparse(truss.Bars(:,1),truss.Bars(:,2),abs(IntensityDataM(:,k)),size(truss.Node,1),size(truss.Node,1));
    VIntensityDataInten(:,k) = sum((IntensityDataIntenk+IntensityDataIntenk'),2); 
end
%VisualFold(Uhis(:,1:interv:endicrm),truss,angles,Fhis(1:interv:endicrm,:),instdof,'IntensityMap','Vertex','IntensityData',VIntensityDataInten)
% Animation monitoring panel-wise value change
%VisualFold(Uhis(:,1:10:endicrm),truss,angles,[],[],'IntensityMap','Edge','IntensityData',STAT.bar.Sx(:,1:10:endicrm),'ShowInitial','off')
% Animation only showing the configurational change
VisualFold(Uhis(:,1:interv:endicrm),truss,angles,[],[])

%% Inspect nodal index assignment
figure()
PlotOri(Node,Panel,[],'ShowNumber','on');
axis equal

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
    dsp = sign(instdof(2))*Uhis((instdof(1)*3-(3-abs(instdof(2)))),:);
    subplot(2,2,4);
    dsp = sign(instdof(2))*Uhis((instdof(1)*3-(3-abs(instdof(2)))),:);
    plot(dsp(1:IraIcrm),Fhis(1:IraIcrm),'-','Color','k','LineWidth',1.5)
    hold on
    plot(dsp(IraIcrm:end),Fhis(IraIcrm:end),'--','Color','k','LineWidth',1.5)
    axis tight
    xlabel('\eta','fontsize',12)
    ylabel('\epsilon','fontsize',12); 
    
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

% Stored energy vs Displacement (advanced)
figure()
plot(dsp,STAT.PE,'r-','linewidth',2);    % Red line is the total energy.
hold on                                  % Between red and cyan is the folding energy. 
plot(dsp,STAT.bend.UB+STAT.bar.US,'c-'); % Between cyan and magenta is the portion of energy for bending.
plot(dsp,STAT.bar.US,'m-');              % Below magenta is the stretching energy of bars.
axis tight
xlabel('Displacement','fontsize',14);
ylabel('Stored Energy','fontsize',14);

%% Export final configuration to an OBJ file
% Write2OBJ('BendedMiura5x5', Nodew, truss.Trigl, truss.Bars, angles.bend)