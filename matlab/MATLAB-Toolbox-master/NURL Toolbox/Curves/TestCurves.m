clear; clc;

lin1=nrlline([0,0], [1,0]);
lin2=nrlline([0,0.6], [1,1.1]);

nrlplot(lin1);
hold on;
nrlplot(lin2);

[crvs1, center1, crvs2, center2, crvs]=tangentarcs(lin1, lin2);

figure; hold on;
for i=1:3
    nrlplot(crvs1(i));
end
plot(center1(1), center1(2), 'ro');

figure; hold on;
for i=1:3
    nrlplot(crvs2(i));
end
plot(center2(1), center2(2), 'ro');

figure; hold on;
for i=1:4
    nrlplot(crvs(i));
end
plot(center1(1), center1(2), 'ro');
plot(center2(1), center2(2), 'ro');







