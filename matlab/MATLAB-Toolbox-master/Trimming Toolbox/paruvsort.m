function [paru,parv]=paruvsort(u,v)
% sort the order of parameters u,v in the corresponding direction
% input:
%       u:parameter coordinates when v in y-direction keeps constant of two surfaces(uv1)
%       v:parameter coordinates when u in x-direction keeps constant of two surfaces(uv2)
%output:
%       paru:paramater coordinates from left to right in x-direction in
%       parameter domain of surface 1(sphere)
%       parv:paramater coordinates from left to right in x-direction according to the order of surface 1 in
%       parameter domain of surface 2(quadratic surface)

A=v{1,:};
B=v{2,:};

[~,n]=size(u);
for i=1:n
    if (~isempty(u{1,i}))
        A=[A,u{1,i}];
    end
     if (~isempty(u{2,i}))
        B=[B,u{2,i}];
     end
end

[~,j]=sort(A(1,:),2);
paru=A(:,j);
parv=B(:,j);
end