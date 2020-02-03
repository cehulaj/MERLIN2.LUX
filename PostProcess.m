function STAT = PostProcess(U_his,truss,angles,IraIcrm)
%% Get Data
Exbar = zeros(size(truss.Bars,1),size(U_his,2)); 
FdAngle = zeros(size(angles.fold,1),size(U_his,2)); 
BdAngle = zeros(size(angles.bend,1),size(U_his,2));

eDofb = kron(truss.Bars,3*ones(1,3))+repmat([-2,-1,0],size(truss.Bars,1),2);
for icrm=1:size(U_his,2)
    Ui = U_his(:,icrm);
    Nodenw = truss.Node;
    Nodenw(:,1) = truss.Node(:,1)+Ui(1:3:end);
    Nodenw(:,2) = truss.Node(:,2)+Ui(2:3:end);
    Nodenw(:,3) = truss.Node(:,3)+Ui(3:3:end);
    
    du = [Ui(eDofb(:,1:3))-Ui(eDofb(:,4:6))];
    Exbar(:,icrm) = truss.B*Ui./truss.L+0.5*sum(du.^2,2)./(truss.L.^2);

    for del = 1:size(angles.bend,1)
        bend = angles.bend(del,:);
        BdAngle(del,icrm) = FoldKe(Nodenw,bend);
    end

    for fel = 1:size(angles.fold,1)
        fold = angles.fold(fel,:);
        FdAngle(fel,icrm) = FoldKe(Nodenw,fold);
    end
end

% Modified
Exbar2 = zeros(size(truss.Bars,1),size(U_his,2));
UDef = U_his(:,IraIcrm);
for icrm2 = IraIcrm:size(U_his,2)
    Ui2 = U_his(:,icrm2) - UDef;
    
    du2 = [Ui2(eDofb(:,1:3))-Ui2(eDofb(:,4:6))];
    Exbar2(:,icrm2) = truss.B2*Ui2./truss.L2+0.5*sum(du2.^2,2)./(truss.L2.^2);
end    

%% Interpret Data
% Modified
[Sx_bar, ~, Wb] = truss.CM(Exbar);
[Sx_bar2, ~, Wb2] = truss.CM2(Exbar2);

Rspr_fd = zeros(size(FdAngle)); Efold = Rspr_fd;
Rspr_bd = zeros(size(BdAngle)); Ebend = Rspr_bd;

% Modified
for i = 1:(IraIcrm-1)
    [Rspr_fdi, ~, Efoldi] = angles.CMfold(FdAngle(:,i),angles.pf0,angles.Kf,truss.L((size(angles.bend,1)+1):(size(angles.bend,1)+size(angles.fold,1))));
    [Rspr_bdi, ~, Ebendi] = angles.CMbend(BdAngle(:,i),angles.pb0,angles.Kb,truss.L(1:size(angles.bend,1)));
    Rspr_fd(:,i) = Rspr_fdi; Efold(:,i) = Efoldi;
    Rspr_bd(:,i) = Rspr_bdi; Ebend(:,i) = Ebendi;
end

% Modified
for i = IraIcrm:size(U_his,2)
    [Rspr_fdi, ~, Efoldi] = angles.CMfold(FdAngle(:,i),angles.pf0,angles.Kf,truss.L((size(angles.bend,1)+1):(size(angles.bend,1)+size(angles.fold,1))));
    [Rspr_bdi, ~, Ebendi] = angles.CMbend(BdAngle(:,i),angles.pb0,angles.Kb,truss.L(1:size(angles.bend,1)));
    [Rspr_fdi2, ~, Efoldi2] = angles.CMfold2(FdAngle(:,i),angles.pf02,angles.Kf2,truss.L((size(angles.bend,1)+1):(size(angles.bend,1)+size(angles.fold,1))));
    [Rspr_bdi2, ~, Ebendi2] = angles.CMbend2(BdAngle(:,i),angles.pb02,angles.Kb2,truss.L(1:size(angles.bend,1)));
    Rspr_fd(:,i) = Rspr_fdi+Rspr_fdi2; Efold(:,i) = Efoldi+Efoldi2;
    Rspr_bd(:,i) = Rspr_bdi+Rspr_bdi2; Ebend(:,i) = Ebendi+Ebendi2;
end

% Modified
STAT.bar.Ex1 = Exbar;
STAT.bar.Ex2 = Exbar2;
STAT.bar.Sx = Sx_bar+Sx_bar2;    
STAT.bar.USi = diag(truss.L.*truss.A)*(Wb+Wb2);
STAT.bar.US = sum(STAT.bar.USi,1);

STAT.fold.Angle = FdAngle;
STAT.fold.RM = Rspr_fd;
STAT.fold.UFi = Efold;
STAT.fold.UF = sum(STAT.fold.UFi,1);

STAT.bend.Angle = BdAngle;
STAT.bend.RM = Rspr_bd;
STAT.bend.UBi = Ebend;
STAT.bend.UB = sum(STAT.bend.UBi,1);

STAT.PE = STAT.bar.US+STAT.fold.UF+STAT.bend.UB;

