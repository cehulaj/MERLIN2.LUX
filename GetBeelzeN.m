function [Node, Panel] = GetBeelzeN(N,R1,R2,R3)
Node = zeros(3*N+1,3);
for n = 1:N
    Node(n,1) = (R1+R2)*cos((n-1)*2*pi/N);
    Node(n,2) = (R1+R2)*sin((n-1)*2*pi/N);
    Node(N+n,1) = R1*cos((n-1)*2*pi/N);
    Node(N+n,2) = R1*sin((n-1)*2*pi/N);
    Node(2*N+n,1) = R3*cos((n-1/2)*2*pi/N);
    Node(2*N+n,2) = R3*sin((n-1/2)*2*pi/N);
end

PMat = zeros(4*N,3);
for n = 1:N
    PMat(n,1) = n;
    PMat(n,2) = 2*N+n;
    PMat(n,3) = N+n;
    PMat(N+n,1) = N+n;
    PMat(N+n,2) = 2*N+n;
    PMat(N+n,3) = 3*N+1;
    PMat(2*N+n,1) = n+1;
    PMat(2*N+n,2) = N+n+1;
    PMat(2*N+n,3) = 2*N+n;
    PMat(3*N+n,1) = N+n+1;
    PMat(3*N+n,2) = 3*N+1;
    PMat(3*N+n,3) = 2*N+n;
end
PMat(3*N,1) = 1;
PMat(3*N,2) = N+1;
PMat(3*N,3) = 3*N;
PMat(4*N,1) = N+1;
PMat(4*N,2) = 3*N+1;
PMat(4*N,3) = 3*N;
Panel = mat2cell(PMat,ones(size(PMat,1),1),3);

end