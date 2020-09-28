%% Compute the nonvanishing basis functions & their derivatives
% Input:  u, knotspanIndex, degree, knotVector
% Output: the basis function values & their derivatives
% only works for knotVector has no double knots or more except first & last 3 terms
% The number of knots is always equal to the number of control points plus curve degree plus one (i.e. number of control points plus curve order).
% the half-open interval, [u;, ui+1 ), special case when u = ui+1.
% Xu Yi, 2019

%%
function [basisFunValue, dersBasisFunValue] = BasisFuns(u, degree, knotVector)
% basisFunValue N(u) % knotspanIndex = i+1;
knotspanIndex = FindSpan(u, knotVector); % returns the span index.
% computes all the nonzero basis functions @u and stores them in the array N[i,:].
N = zeros(degree+1,degree+1); N(1,degree+1) = 1; % initialize. Ndegree+1,0 = 1
for i = 1:degree % order = degree+1;
    for j = 1:i+1 % the NURBS book 2nd. P70
        left1 = u - knotVector( knotspanIndex-( i-j ) );
        left2 = u - knotVector( knotspanIndex-( i+1-j ) );
        right1 = knotVector( knotspanIndex-1+j ) - u;
        if knotspanIndex+j > length(knotVector) % special case
            right2 = 0;
        else
            right2 = knotVector( knotspanIndex+j ) - u;
        end
        Para1 = left2/(right1+left2);
        Para2 = right2/(right2+left1);
        if right1+left2 == 0 % define 0/0=0
            Para1 = 0;
        elseif right2+left1 == 0
            Para2 = 0;
        end
        Nindex = degree+1-(i+1-j); % N: Lower right triangle
        if j == 1 % the first term of Eq. (2.14) are not computed, since Ni = 0;
            N(i+1,Nindex) = Para2 * N(i,Nindex+1);
        elseif j == i+1 % the last term of Eq. (2.16) are not computed, since Ni = 0;
            N(i+1,Nindex) = Para1 * N(i,Nindex);
        else
            N(i+1,Nindex) = Para1 * N(i,Nindex) + Para2 * N(i,Nindex+1);
        end
    end
end
% derivatives % the NURBS book 2nd. P61 (2.10)
% parameter a
Para_a = zeros(degree+1,degree+1,degree+1); Para_a(1,1,:) = 1; % initialize. a0,0 = 1
for i = 1:degree % a0,0 is initialized / ith derivative
    for j = 1:i+1
        for k = 1:degree+1 % a is subject to i(which in Ni)
            if knotspanIndex+k-(degree+1) + degree +j-1 -i+1 > length(knotVector) % special case
                temp = 0;
            else
                temp = knotVector( knotspanIndex+k-(degree+1) + degree +j-1 -i+1 ) - knotVector( knotspanIndex+k-(degree+1) +j-1 );
            end
            if j == 1
                Para_a(i+1,j,k) = Para_a(i,1,k) / temp;
            elseif j == i+1
                Para_a(i+1,j,k) = - Para_a(i,i,k) / temp;
            else
                Para_a(i+1,j,k) = ( Para_a(i,j,k) - Para_a(i,j-1,k) ) / temp;
            end
            if temp == 0 % special case. a/0 = 0
                Para_a(i+1,j,k) = 0;
            end
        end
    end
end
% dersBasisFunValue / calculate derivatives of Ni,degree
ders = zeros(degree,degree+1); % initialize. [der x degree+1]
for i = 1:degree % derivative. k<=p
    for j = 1:degree+1
        temp1 = factorial(degree) / factorial(degree-i);
        temp2 = 0;
        for k = 1:i+1
            temp211 = j+k-1;
            if temp211 > degree+1
                temp21 = 0;
            elseif temp211 < 1
                temp21 = 0;
            else
                temp21 = N( degree+1-i,temp211 );
            end
            temp2 = temp2 + Para_a(i+1,k,j)*temp21 ;
        end
        ders(i,j) = temp1 * temp2;
    end
end

basisFunValue = N;
dersBasisFunValue = ders;
end
