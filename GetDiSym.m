function [Node, Panel] = GetDiSym(N,R,h,lyr,phi)
if numel(phi) == 1
    rotangle = (0:lyr-1)*phi;
else
    rotangle = phi;
end
rdl = zeros(lyr,N);
for i = 1:lyr
    rdl(i,:) = linspace(rotangle(i),2*pi/N*(N-1)+rotangle(i),N);
end
Xcood = R*cos(reshape(rdl',[],1));
Ycood = R*sin(reshape(rdl',[],1));
Zcood = h*reshape(repmat(0:lyr-1,N,1),[],1);
Node = [Xcood,Ycood,Zcood];

PMat = zeros(2*N*(lyr-1),3);
for i = 1:lyr-1
    PMat((1:N)+2*(i-1)*N,:) = [1:N;(1:N)+N;mod(1:N,N)+N+1]'+(i-1)*N;
    PMat((N+1:2*N)+2*(i-1)*N,:) = [1:N;mod(1:N,N)+1;mod(1:N,N)+N+1]'+(i-1)*N;
end
Panel = mat2cell(PMat,ones(size(PMat,1),1),3);
