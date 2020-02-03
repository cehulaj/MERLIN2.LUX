function [Uhis,Fhis,angles,truss,IraIcrm] = PathAnalysis(truss,angles,AnalyInputOpt)
tol = 1e-6; MaxIter = 50; 
Node = truss.Node;
if ~isfield(truss,'U0'), truss.U0 = zeros(3*size(truss.Node,1),1); end  
U = truss.U0;

if strcmpi(AnalyInputOpt.LoadType, 'Force')
    MaxIcr = AnalyInputOpt.MaxIcr;
    b_lambda = AnalyInputOpt.InitialLoadFactor;
    Uhis = zeros(3*size(Node,1),MaxIcr);
    FreeDofs = setdiff(1:3*size(Node,1),truss.FixedDofs);
    lmd = 0; icrm = 0; MUL = [U,U];
    Fhis = zeros(MaxIcr,1);
    if isfield(AnalyInputOpt,'Load')
        F = AnalyInputOpt.Load;
    end
    
    % Modified
    while lmd<1 && icrm<MaxIcr && ~AnalyInputOpt.StopCriterion(Node,U,icrm,lmd)
        icrm = icrm+1;
        iter = 0; err = 1;
        %fprintf('icrm = %d, lambda = %6.4f\n',icrm,lmd);
        
        %Modified
        if isfield(AnalyInputOpt,'AdaptiveLoad')
            F = AnalyInputOpt.AdaptiveLoad(Node,U,icrm,AnalyInputOpt);
        end
        while err>tol && iter<MaxIter
            iter = iter+1;
            [IF,K] = GlobalK_fast_ver(U,Node,truss,angles);
            R = lmd*F-IF;   MRS = [F,R];
            MUL(FreeDofs,:) = K(FreeDofs,FreeDofs)\MRS(FreeDofs,:);
            dUp = MUL(:,1); dUr = MUL(:,2);
            if iter==1, dUr = 0*dUr; end
            dlmd=nlsmgd(icrm,iter,dUp,dUr,b_lambda);
            dUt = dlmd*dUp+dUr;
            U = U+dUt;
            err = norm(dUt(FreeDofs));
            lmd = lmd+dlmd;
            %fprintf('    iter = %d, err = %6.4f, dlambda = %6.4f\n',iter,err,dlmd);
            if err > 1e8, disp('Divergence!'); break; end
        end

        if iter>15
            b_lambda = b_lambda/2;
            disp('Reduce constraint radius...')
            icrm = icrm-1;
            U = Uhis(:,max(icrm,1));  % restore displacement
            lmd = Fhis(max(icrm,1));   % restore load
        elseif iter<3
            disp('Increase constraint radius...')
            b_lambda = b_lambda*1.5;
            Uhis(:,icrm) = U;
            Fhis(icrm) = lmd; 
        else
            Uhis(:,icrm) = U;
            Fhis(icrm) = lmd; 
        end
    end
    
    fprintf('icrm = %d, lambda = %6.4f\n',icrm,lmd);
    
    % Modified
    % Pause for irradiation
    Pause = round(MaxIcr/10);
    for i = 1:Pause
        icrm = icrm+1;
        Uhis(:,icrm) = U;
        Fhis(icrm) = lmd;
    end
    % When did the irradiation happened?
    IraIcrm = icrm;
    % Lmd will decrease
    b_lambda = -b_lambda;
    % Displacement during the irradiation 
    UDef = U;
    % New reference configuration
    NodeDef(:,1) = Node(:,1)+U(1:3:end);
    NodeDef(:,2) = Node(:,2)+U(2:3:end);
    NodeDef(:,3) = Node(:,3)+U(3:3:end);
    % New neutral angles    
    pf02 = zeros(size(angles.fold,1),1);        
    for i = 1:size(angles.fold,1)
        pf02(i) = FoldKe(NodeDef,angles.fold(i,:));
    end
    pb02 = zeros(size(angles.bend,1),1);    
    for i = 1:size(angles.bend,1)
        pb02(i) = FoldKe(NodeDef,angles.bend(i,:));
    end    
    angles.pf02 = pf02;
    angles.pb02 = pb02;
    % New compatibility matrix and lengths of bar elements
    [B2, L2] = dirc3d(NodeDef,truss.Bars);    
    truss.B2 = B2; 
    truss.L2 = L2;       
    
    % Modified 
    % Unload
    while lmd>0 && icrm<MaxIcr && ~AnalyInputOpt.StopCriterion(Node,U,icrm,lmd)
        icrm = icrm+1;
        iter = 0; err = 1;
        %fprintf('icrm = %d, lambda = %6.4f\n',icrm,lmd);
       
        % Modified
        if isfield(AnalyInputOpt,'AdaptiveLoad')
            F = AnalyInputOpt.AdaptiveLoad(Node,U,icrm,AnalyInputOpt);
        end
        while err>tol && iter<MaxIter
            iter = iter+1;
            
            % Modified
            [IF,K] = GlobalK_fast_ver2(U,UDef,Node,truss,angles);
            
            R = lmd*F-IF;   MRS = [F,R];
            MUL(FreeDofs,:) = K(FreeDofs,FreeDofs)\MRS(FreeDofs,:);
            dUp = MUL(:,1); dUr = MUL(:,2);
            if iter==1, dUr = 0*dUr; end
            dlmd=nlsmgd(icrm,iter,dUp,dUr,b_lambda);
            dUt = dlmd*dUp+dUr;
            U = U+dUt;
            err = norm(dUt(FreeDofs));
            lmd = lmd+dlmd;
            %fprintf('    iter = %d, err = %6.4f, dlambda = %6.4f\n',iter,err,dlmd);
            if err > 1e8, disp('Divergence!'); break; end
        end

        if iter>15
            b_lambda = b_lambda/2;
            disp('Reduce constraint radius...')
            icrm = icrm-1;
            U = Uhis(:,max(icrm,1));  % restore displacement
            lmd = Fhis(max(icrm,1));   % restore load
        elseif iter<3
            disp('Increase constraint radius...')
            b_lambda = b_lambda*1.5;
            Uhis(:,icrm) = U;
            Fhis(icrm) = lmd; 
        else
            Uhis(:,icrm) = U;
            Fhis(icrm) = lmd; 
        end
    end
    fprintf('icrm = %d, lambda = %6.4f\n',icrm,lmd);
    % Modified
    % Pause for the video
    Pause = round(MaxIcr/20);
    for i = 1:Pause
        icrm = icrm+1;
        Uhis(:,icrm) = U;
        Fhis(icrm) = lmd;
    end
    
elseif strcmpi(AnalyInputOpt.LoadType, 'Displacement')
    % Modified
    MaxIcr = AnalyInputOpt.MaxIcr;
    % Rewritten: Uhis = zeros(3*size(Node,1),AnalyInputOpt.DispStep*2);
    Uhis = zeros(3*size(Node,1),MaxIcr);
    
    if isfield(AnalyInputOpt,'Load')
        Fdsp = AnalyInputOpt.Load/AnalyInputOpt.DispStep;
    else
        % Modified
        Fdsp = AnalyInputOpt.AdaptiveLoad(Node,U,1,AnalyInputOpt);
    end
    ImpDofs = find(Fdsp~=0);
    FreeDofs = setdiff(setdiff(1:3*size(Node,1),truss.FixedDofs),ImpDofs);
    icrm = 0;  dspmvd = 0;  attmpts = 0;
    mvstepsize = 1;  damping = 1;
    % Rewritten: Fhis = zeros(AnalyInputOpt.DispStep,numel(ImpDofs));
    Fhis = zeros(MaxIcr,numel(ImpDofs));
    while (dspmvd <= 1 && ~AnalyInputOpt.StopCriterion(Node,U,icrm,0.5)) && attmpts <= 20
        icrm = icrm+1;
        iter = 0; err = 1;   
        %fprintf('icrm = %d, dspimps = %6.4f\n',icrm,dspmvd);
        % Modified
        if isfield(AnalyInputOpt,'AdaptiveLoad')
            Fdsp = AnalyInputOpt.AdaptiveLoad(Node,U,icrm,AnalyInputOpt);
        end
        U = U+mvstepsize*Fdsp;
        U(truss.FixedDofs)=0;
        while err>tol && iter<MaxIter
            iter = iter+1;
            [IF,K] = GlobalK_fast_ver(U,Node,truss,angles);
            dU = zeros(3*size(Node,1),1);
            dU(FreeDofs) = K(FreeDofs,FreeDofs)\(-IF(FreeDofs));
            err = norm(dU(FreeDofs));
            U = U+damping*dU; 
            %fprintf('    iter = %d, err = %6.4f\n',iter,err);
        end

        if iter>=((mvstepsize>1)+1)*MaxIter/(damping+1)  
            % an aggressive step needs more iterations
            attmpts = attmpts+1;
            icrm = icrm-1;
            if attmpts<=10  
                mvstepsize = mvstepsize*0.5; 
                disp('Take a more conservative step...')
            else
                mvstepsize = max(mvstepsize,1)*1.5;  
                damping = damping*0.75;
                disp('Take a more aggressive step...')
            end
            U = Uhis(:,max(icrm,1)); % restore displacement            
        else
            dspmvd = dspmvd+mvstepsize/AnalyInputOpt.DispStep;
            attmpts = 0;
            damping = 1;
            if mvstepsize<1
                mvstepsize = min(mvstepsize*1.1,1); % gradually go back to 1
            else
                mvstepsize = max(mvstepsize*0.9,1);
            end
            Uhis(:,icrm) = U;
            [Fend,~] = GlobalK_fast_ver(U,Node,truss,angles);
            Fhis(icrm,:) = -Fend(ImpDofs)'; 
        end
    end
    
    fprintf('icrm = %d, dspimps = %6.4f\n',icrm,dspmvd);
    
    % Modified
    FreeDofs = setdiff(1:3*size(Node,1),truss.FixedDofs);
    if size(Fhis,2)>1, Fhis = sum(Fhis,2); end
    lmd = Fhis(icrm);
    F = Fend/lmd;
    MUL = [U,U];
    % Pause for irradiation
    Pause = round(MaxIcr/10);
    for i = 1:Pause
        icrm = icrm+1;
        Uhis(:,icrm) = U;
        Fhis(icrm) = lmd;
    end
    % When did the irradiation happened?
    IraIcrm = icrm;
    % Lmd will approach zero
    b_lambda = -sign(lmd)*AnalyInputOpt.InitialLoadFactor;
    % Displacement during the irradiation 
    UDef = U;
    % New reference configuration
    NodeDef(:,1) = Node(:,1)+U(1:3:end);
    NodeDef(:,2) = Node(:,2)+U(2:3:end);
    NodeDef(:,3) = Node(:,3)+U(3:3:end);
    % New neutral angles    
    pf02 = zeros(size(angles.fold,1),1);        
    for i = 1:size(angles.fold,1)
        pf02(i) = FoldKe(NodeDef,angles.fold(i,:));
    end
    pb02 = zeros(size(angles.bend,1),1);    
    for i = 1:size(angles.bend,1)
        pb02(i) = FoldKe(NodeDef,angles.bend(i,:));
    end    
    angles.pf02 = pf02;
    angles.pb02 = pb02;
    % New compatibility matrix and lengths of bar elements
    [B2, L2] = dirc3d(NodeDef,truss.Bars);    
    truss.B2 = B2; 
    truss.L2 = L2;
    
    % Modified
    while sign(Fhis(IraIcrm))*lmd>0 && icrm<MaxIcr && ~AnalyInputOpt.StopCriterion(Node,U,icrm,lmd)
        icrm = icrm+1;
        iter = 0; err = 1;
        %fprintf('icrm = %d, lambda = %6.4f\n',icrm,lmd);
        
        % Modified
        if isfield(AnalyInputOpt,'AdaptiveLoad')
            F = AnalyInputOpt.AdaptiveLoad(Node,U,icrm,AnalyInputOpt);
        end
        while err>tol && iter<MaxIter
            iter = iter+1;
            
            % Modified
            [IF,K] = GlobalK_fast_ver2(U,UDef,Node,truss,angles);
            
            R = lmd*F-IF;   MRS = [F,R];
            MUL(FreeDofs,:) = K(FreeDofs,FreeDofs)\MRS(FreeDofs,:);
            dUp = MUL(:,1); dUr = MUL(:,2);
            if iter==1, dUr = 0*dUr; end
            dlmd=nlsmgd(icrm-IraIcrm,iter,dUp,dUr,b_lambda);
            dUt = dlmd*dUp+dUr;
            U = U+dUt;
            err = norm(dUt(FreeDofs));
            lmd = lmd+dlmd;
            %fprintf('    iter = %d, err = %6.4f, dlambda = %6.4f\n',iter,err,dlmd);
            if err > 1e8, disp('Divergence!'); break; end
        end

        if iter>15
            b_lambda = b_lambda/2;
            disp('Reduce constraint radius...')
            icrm = icrm-1;
            U = Uhis(:,max(icrm,1));  % restore displacement
            lmd = Fhis(max(icrm,1));   % restore load
        elseif iter<3
            disp('Increase constraint radius...')
            b_lambda = b_lambda*1.5;
            Uhis(:,icrm) = U;
            Fhis(icrm) = lmd; 
        else
            Uhis(:,icrm) = U;
            Fhis(icrm) = lmd; 
        end
    end
    fprintf('icrm = %d, lambda = %6.4f\n',icrm,lmd);
    % Modified
    % Pause for the video
    Pause = round(MaxIcr/20);
    for i = 1:Pause
        icrm = icrm+1;
        Uhis(:,icrm) = U;
        Fhis(icrm) = lmd;
    end
    
else
    disp('Unknown load type!!!')
end

icrm = icrm+1;
Uhis(:,icrm:end) = [];
Fhis(icrm:end,:) = [];
end

%--------------------------------------------------------------------------
function dl=nlsmgd(step,ite,dup,dur,cmp)
% Modified Generalized Displacement Control Method
global dupp1 sinal dupc1 numgsp
if ite==1
    if step==1
        sinal=sign(dot(dup,dup));
        dl=cmp;
        numgsp=dot(dup,dup);   
    else
        sinal=sinal*sign(dot(dupp1,dup));
        gsp=numgsp/dot(dup,dup);
        dl=sinal*cmp*sqrt(gsp);
    end 
    dupp1=dup;
    dupc1=dup;
else
    dl=-dot(dupc1,dur)/dot(dupc1,dup);
end
end

%--------------------------------------------------------------------------
function [B, L] = dirc3d(Node,Ele)
Ne = size(Ele,1); Nn = size(Node,1);
D = [Node(Ele(:,2),1)-Node(Ele(:,1),1), Node(Ele(:,2),2)-Node(Ele(:,1),2),...
    Node(Ele(:,2),3)-Node(Ele(:,1),3)];
L = sqrt(D(:,1).^2+D(:,2).^2+D(:,3).^2);
D = [D(:,1)./L D(:,2)./L D(:,3)./L];
B = sparse(repmat((1:Ne)',1,6),[3*Ele(:,1)-2 3*Ele(:,1)-1 3*Ele(:,1),...
           3*Ele(:,2)-2 3*Ele(:,2)-1 3*Ele(:,2)],[D -D],Ne,3*Nn);
B = -B;
end

