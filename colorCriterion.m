function a = colorCriterion(M,N)
% Inner product of two images M,N
%Copyright: Wu Shiqian
%8 May 2009
if isa(M,'uint8')
    M = double(M);
end
if isa(N,'uint8')
    N = double(N);
end
M = log(M+1);
N = log(N+1);
a = sqrt(sum((M+eps).*(M+eps),3)).*sqrt(sum((N+eps).*(N+eps),3));
a = dot(M+eps,N+eps,3)./a;
