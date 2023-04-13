function B = probabilityValue(A)
% A is a vector
% Output B is the probability value from A
%Copyright: Wu Shiqian
%15 May 2009
if isa(A,'uint8')
    A = double(A);
end
if isempty (A)
    disp('A is a empty vector')
    return
end
sig =1;
T = length(A);
B = [];
while sig
    tmp = A(1);
    n = find(A==tmp);
    b1 = [tmp,length(n)];
    B = [B;b1];
    A(n)=[];
    if isempty(A)
        sig =0;
    end
end
B = sum(B(:,1).*B(:,2))/T;