

Q=Qxy;
[m, n]=size(Q);
for i=1:m
    for j=1:n
        if abs(Q(i,j))<1e-12
            Q(i,j)=0; 
        end
    end
end
sQ=sparse(Q);
figure; spy(sQ); 