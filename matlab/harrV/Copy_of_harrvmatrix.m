clear
k=3;
n=3;
x=1/(2^(n+1)*(k+1)):1/(2^n*(k+1)):((2^(n+1)*(k+1))-1)/(2^(n+1)*(k+1));
x0=1/(2*(k+1)):1/(k+1):(2*(k+1)-1)/(2*(k+1));
y{1}=(2^n)^0.5*3^0.5*(1-2*x0);
y{2}=(2^n)^0.5*5^0.5*(6*x0.^2-6*x0+1);
y{3}=(2^n)^0.5*7^0.5*(-20*x0.^3+30*x0.^2-12*x0+1);

for i=1:n+1
%     P{i}=hadamard(2^(i-1));
    T{i}=kron(eye(2^n/2^(i-1)),hadamard(2^(i-1)));
end
% for j=1:n+1
%     T{j}=P{j};
% for i=1:2^n/2^(j-1)-1
%     T{j}=blkdiag(T{j},P{j});
% end
% end

% j=2;
for j=1:k
    for l=1:2^n
        for i=1:2^n
            V(l+(j-1)*2^n,(1+(i-1)*(k+1)):(k+1+(i-1)*(k+1)))=y{j}*T{1}(l,i);
        end
    end
end
% A=y{1}*T{1}(1,2);
% for i=1:size(V,1)
% % NP=normalize(A,'norm',2);
%     VN(i,:)=normalize(V(i,:),'norm',2);
% end
% plot(x,V(17,:))
% save matrixname VN
har=[ones(1,k+1),-ones(1,k+1)];
H{1}=ones(1,2^n*(k+1));
for i=2:n+1
    H{i}=kron(kron(eye(2^(i-2)),har),ones(1,2^(n-i+1)));
end
Harr=cell2mat(H');
HarrV=[Harr;V];
for i=1:size(HarrV,1)
% NP=normalize(A,'norm',2);
    HarrVN(i,:)=normalize(HarrV(i,:),'norm',2);
end
