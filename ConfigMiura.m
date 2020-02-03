function [NODE,PANEL,BDRY] = ConfigMiura(x_divs,y_divs,theta,a,b,gmma)
theta=pi*theta/180;
gmma = gmma*pi/180;
numx = 2*x_divs; numy = 2*y_divs;
h = a*sin(theta)*sin(gmma);
s = b*cos(gmma)*tan(theta)/sqrt(1+cos(gmma)^2*tan(theta)^2);
l = a*sqrt(1-sin(gmma)^2*sin(theta)^2);
v = b/sqrt(1+cos(gmma)^2*tan(theta)^2);
X = repmat(s*(0:1:numx),numy+1,1);
Y = repmat(l*[0:1:numy]',1,numx+1);
Y(:,2:2:end) = Y(:,2:2:end)+v;
Z = 0*X;
Z(2:2:end,:) = h;
NODE = [reshape(X,[],1),reshape(Y,[],1),reshape(Z,[],1)];

if nargout > 2
    BDRY = zeros(2*numx+2*numy,2);
    Node_idx = reshape(1:size(NODE,1),size(X));
    Lbdry = Node_idx(:,1);
    Rbdry = Node_idx(:,end);
    Bbdry = Node_idx(1,:);
    Tbdry = Node_idx(end,:);

    % Specify boundary
    BDRY(1:numel(Lbdry)-1,:) = [Lbdry(1:end-1) Lbdry(2:end)];
    count = numel(Lbdry)-1;
    BDRY(count+1:count+numel(Bbdry)-1,:) = [Bbdry(1:end-1)' Bbdry(2:end)'];
    count = count+numel(Bbdry)-1;
    BDRY(count+1:count+numel(Rbdry)-1,:) = [Rbdry(1:end-1) Rbdry(2:end)];
    count = count+numel(Rbdry)-1;
    BDRY(count+1:count+numel(Tbdry)-1,:) = [Tbdry(1:end-1)' Tbdry(2:end)'];
    count = count+numel(Tbdry)-1;

    k = 0; PANEL = cell(numx*numy,1);
    for j=1:numy, for i=1:numx
            k = k+1;
            n1 = (i-1)*(numy+1)+j; n2 = i*(numy+1)+j;
            PANEL{k} = [n1 n2 n2+1 n1+1];
    end, end
end

