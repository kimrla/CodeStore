function [crvs, center]=tangentarc(lin1, lin2, m)

%  Get the tangent arcs of two straight lines
% 
%   Inputs: 
%
%      crv1, crv2 - the two curves to be used
% 
%       m - the small (1) or large (2) arc
% 
%  Output: 
%
%       crvs - the resulted curves
% 
%       center - center of the tangent arc
% 
%  Examples:  
% 
%     lin1=nrlline([0,0], [1,0]);
%     lin2=nrlline([0,0.6], [1,1.1]);
%     [crvs1, center1]=tangentarc(lin1, lin2, 1);
%     [crvs2, center2]=tangentarc(lin1, lin2, 2);   
%

% Get the curvature of the lines to judge whether they are straight
t=linspace(0, 1, 20);
[~, cvt1]=curvature(lin1, t); 
[~, cvt2]=curvature(lin2, t); 
c1=max(cvt1); c2=max(cvt2); 
dist1 = nrlmeasure(lin1); 
dist2 = nrlmeasure(lin2); 
if max(abs([c1/dist1, c2/dist2]))>1e-3
    error('One or both of the two lines are not straight.');
end

% Extract the ends and the tangent
pnts1=nrlcrvextract(lin1);
pnts2=nrlcrvextract(lin2);

% Creat analytical straight lines and get the intersection point
ABC1= planeline('2points', pnts1(1:2, 1)', pnts1(1:2, 2)');
ABC2= planeline('2points', pnts2(1:2, 1)', pnts2(1:2, 2)');
inter=planeline('intersection', ABC1, ABC2 );

% Rearrange the direction of the lines
dm1 = DistanceMatrix([inter, 0], pnts1'); 
dm2 = DistanceMatrix([inter, 0], pnts2'); 
[~, i]=min(dm1);
[~, j]=min(dm2);
if i==2
    lin1 = nrlreverse(lin1);
end
if j==2
    lin2 = nrlreverse(lin2);
end

% Extract the two end points of the two lines and get 
%     the directions of the two lines at the two points 
[pnts1, directs1]=nrlcrvextract(lin1);
[pnts2, directs2]=nrlcrvextract(lin2);
dr1=directs1(1:2, 1);
dr2=directs2(1:2, 1);

% Get the bisector vector and angle of two plane vectors
[dr, ~, dr1, dr2]=BisectVector(dr1, dr2);
dr12=[dr1; dr2];

switch m
    case 1
        % Get the center and radius of the small arcs
        pnts=[pnts1(:,1), pnts2(:,1)];
        dm = DistanceMatrix([inter, 0], pnts'); 
        [~, i]=max(dm);
        pnt1=pnts(1:2,i)';
        ABC1= planeline('pointnormal', pnt1, dr12(i, :));
        ABC2= planeline('pointtang', inter, dr);
        center=planeline('intersection', ABC1, ABC2 );
        R=sqrt(sum((pnt1-center).^2));

        % Get the other two ponts of the small arc
        if i==1
            dn=dr2; pnt=pnts2(1:2,1)';
        else
            dn=dr1; pnt=pnts1(1:2,1)';
        end
        ABC1= planeline('pointtang', pnt, dn);
        ABC2= planeline('pointnormal', center, dn);
        pnt2=planeline('intersection', ABC1, ABC2 );
        L=sqrt(sum((inter-center).^2));
        pnt3=inter+(L-R)*dr;

        % Ctreat the small arc by the three points
        arc1=nrl3pntsarc(pnt1, pnt3, pnt2);

        % Get the new lines
        if i==1
            lin11=lin1;
            lin21=nrlline(pnt2, pnts2(1:2,2)');
        else 
            lin11=nrlline(pnt2, pnts1(1:2,2)');
            lin21=lin2;
        end
        
        % Get the final curves or tangent circle
        crvs=[lin11, lin21, arc1];
    case 2
        % Get the center and radius of the large arc
        pnts=[pnts1(:,2), pnts2(:,2)];
        dm = DistanceMatrix([inter, 0], pnts'); 
        [~, j]=min(dm);
        pnt1=pnts(1:2,j)';
        ABC1= planeline('pointnormal', pnt1, dr12(j, :));
        ABC2= planeline('pointtang', inter, dr);
        center=planeline('intersection', ABC1, ABC2 );
        R=sqrt(sum((pnt1-center).^2));

        % Get the other two ponts of the arc
        if j==1
            dn=dr2; pnt=pnts2(1:2,2)';
        else
            dn=dr1; pnt=pnts1(1:2,2)';
        end
        ABC1= planeline('pointtang', pnt, dn);
        ABC2= planeline('pointnormal', center, dn);
        pnt2=planeline('intersection', ABC1, ABC2 );
        pnt3=center+R*dr;

        % Ctreat the arc by the three points
        arc2=nrl3pntsarc(pnt1, pnt3, pnt2);

        % Get the new lines
        if j==1
            lin12=lin1;
            lin22=nrlline(pnt2, pnts2(1:2,1)');
        else 
            lin12=nrlline(pnt2, pnts1(1:2,1)');
            lin22=lin2;
        end

        % Get the final curves or tangent circle
        crvs=[lin12, lin22, arc2];    
end

% Remove the curve that is a point
if R/max([dist1, dist2])<1e-3
    crvs(3)=[];
end


%% Demo
% lin1=nrlline([0,0], [1,0]);
% lin2=nrlline([0,0.6], [1,1.1]);
% 
% nrlplot(lin1);
% hold on;
% nrlplot(lin2);
% 
% [crvs1, center1]=tangentarc(lin1, lin2, 1);
% [crvs2, center2]=tangentarc(lin1, lin2, 2);   
% 
% figure; hold on;
% for i=1:3
%     nrlplot(crvs1(i));
% end
% plot(center1(1), center1(2), 'ro');
% 
% figure; hold on;
% for i=1:3
%     nrlplot(crvs2(i));
% end
% plot(center2(1), center2(2), 'ro');



