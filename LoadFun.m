function FRC = LoadFun(Node,U,icrm,AnalyInputOpt)
%% Get the data
FDef = AnalyInputOpt.FDef;
OrigamiType = AnalyInputOpt.OrigamiType;

%% Compute current configuration
Nini = zeros(3*size(Node,1),1);
Nini(1:3:end) = Node(:,1);
Nini(2:3:end) = Node(:,2);
Nini(3:3:end) = Node(:,3); 
    
NDef = Nini + U;

%% Compute the force/displacement vector FRC (F)
if icrm<=0 
    error('Wrong increment!');
else    
    FRC = zeros(3*size(Node,1),1);
    
    if strcmpi(OrigamiType, 'Bird') 
        FRC(1) = - NDef(3);
        FRC(3) = NDef(1);
        FRC = FDef*FRC/sqrt(NDef(1)^2+NDef(3)^2);
    
    elseif strcmpi(OrigamiType, 'Beelze')
        FRC(1) = NDef(3);
        FRC(3) = - NDef(1);
        FRC = FDef*FRC/sqrt(NDef(1)^2+NDef(3)^2);
        
    elseif strcmpi(OrigamiType, 'BeelzeN')
        N = AnalyInputOpt.N;
        K = zeros(3*N,1);
        for i = 1:3*N
            K(i) = NDef(i) - NDef(3*N+i);
        end
        for n = 1:N
            FRC(3*n-2) = -K(3*n)*cos((n-1)*2*pi/N);
            FRC(3*n-1) = -K(3*n)*sin((n-1)*2*pi/N);
            FRC(3*n) = sqrt((K(3*n-2))^2+(K(3*n-1))^2);
        end
        FRC = FDef*FRC/sqrt(K(1)^2+K(2)^2+K(3)^2);
        
    else    
        error('Unknown OrigamiType!');
    end
end    