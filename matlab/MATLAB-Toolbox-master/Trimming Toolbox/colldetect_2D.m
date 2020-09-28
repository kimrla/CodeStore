function logic = colldetect_2D( box1,box2)
% Collision detection of 2 boxes based on the projection theorem (for AABB and OBB) or using the
% intersection detection among all edges and polygons (for NURBS control polygons).
% This function is to determine whether 2 polygons in 2D are intersected.

% Input:
%     box1,box2: Two boxes to be detected, which are represented by the
%       structure obtained by the function boxconstruct.
% Output:
%     logic: The logical scalar meaning 1: if the two boxes are intersected
%       or 0: if the two boxes are NOT intersected.

form1=box1.form;
form2=box2.form;

if strcmp(form1,'Circle') || strcmp(form2,'Circle') % If circle box exists
    if strcmp(form1,form2) % Circle-Circle
        center1=box1.center;
        center2=box2.center;
        r1=box1.radius;
        r2=box2.radius;
        if norm(center1-center2)>r1+r2
            logic=false;
        else
            logic=true;
        end
    else % Circle-Polygons ( including rectangles)
        if strcmp(form1,'Circle')
            center=box1.center;
            r=box1.radius;
            box=box2;
        else
            center=box2.center;
            r=box2.radius;
            box=box1;
        end
        pts=box.points;
        logic=false;
        % The code may be wrong, needed to be modified.
        for i=1:size(pts,1)-1 % The first element is the same as the last one
            [~,~,dist,~]=propt2line(center,pts(i,:),pts(i+1,:));
            if dist<=r
                logic=true;
                break;
            end
        end
    end
    
elseif  strcmp(form1,'POLYGON') || strcmp(form2,'POLYGON')% Polygons exist
    pts1=box1.points; % The first element is the same as the last one
    pts2=box2.points;
    in=inpolygon(pts1(:,1),pts1(:,2),pts2(:,1),pts2(:,2));
    if sum(in)>0
        logic=true;
    else
        in=inpolygon(pts2(:,1),pts2(:,2),pts1(:,1),pts1(:,2));
        if sum(in)>0
            logic=true;
        else
            logic=false;
        end
    end
    % If logic obtained above is false, then it needs to determine further.
    if logic==false
        % Calculate the local coordinates of all vertices of 2 polygons, regarding each separated axis as
        % the new x-axis, and then compare the x-coordinate of the 2 series
        % of vertices in the new co sys (local sys).
        sa1=box1.sa;
        sa2=box2.sa;
        numsa1=length(sa1);
        while true
            for i=1:numsa1
                theta=anglereverse(sa1{i});
                [~,npt1]=cordtrans([pts1,zeros(size(pts1,1),1)],[0 0 0],theta);
                [~,npt2]=cordtrans([pts2,zeros(size(pts2,1),1)],[0 0 0],theta);
                xspan1=[min(npt1(:,1)),max(npt1(:,2))];
                xspan2=[min(npt2(:,2)),max(npt2(:,2))];
                if xspan1(2)<xspan2(1) || xspan2(2)<xspan1(1) % At least one projection indicates NOT intersecting
                    logic=false;
                    break;
                else
                    logic=true;
                end
            end
            if i<numsa1
                break;
            elseif logic==true
                numsa2=length(sa2);
                for i=1:numsa2
                    theta=anglereverse(sa2{i});
                    [~,npt1]=cordtrans([pts1,zeros(size(pts1,1),1)],[0 0 0],theta);
                    [~,npt2]=cordtrans([pts2,zeros(size(pts2,1),1)],[0 0 0],theta);
                    xspan1=[min(npt1(:,1)),max(npt1(:,2))];
                    xspan2=[min(npt2(:,2)),max(npt2(:,2))];
                    if xspan1(2)<xspan2(1) || xspan2(2)<xspan1(1) % At least one projection indicates NOT intersecting
                        logic=false;
                        break;
                    else
                        logic=true;
                    end
                end 
            end
            break;
        end
    end
        
elseif  strcmp(form1,'OBB') || strcmp(form2,'OBB')% OBB exists and Polygons don't exist
    sa1=box1.sa;pts1=box1.points;
    sa2=box2.sa;pts2=box2.points;
    sa={sa1{1},sa1{2},sa2{1},sa2{2}};
    center1=(pts1(1,:)+pts1(3,:))/2;
    center2=(pts2(1,:)+pts2(3,:))/2;
    center=center2-center1;
    h1=(pts1(2,:)-pts1(1,:))/2;
    w1=(pts1(3,:)-pts1(2,:))/2;
    h2=(pts2(2,:)-pts2(1,:))/2;
    w2=(pts2(3,:)-pts2(2,:))/2;
    % Remove repeated sa
    sa=unique(sa,'rows');
    logic=true;
    for i=1:length(sa)
        det1=abs(dot(center,sa{i}));
        det2=max(abs(dot((h1+w1),sa{i})),abs(dot((h1-w1),sa{i})))+...
            max(abs(dot((h2+w2),sa{i})),abs(dot((h2-w2),sa{i})));              
        if det1>det2
            logic=false;
            break;
        end
    end
    
else % AABB-AABB
    logic=true;
    pts1=box1.points;
    pts2=box2.points;
    xspan1=[pts1(1,1),pts1(3,1)];
    xspan2=[pts2(1,1),pts2(3,1)];
    yspan1=[pts1(1,2),pts1(3,2)];
    yspan2=[pts2(1,2),pts2(3,2)];
    if xspan1(2)<xspan2(1) || xspan2(2)<xspan1(1) % At least one projection indicates NOT intersecting
        logic=false;
    elseif yspan1(2)<yspan2(1) || yspan2(2)<yspan1(1) % At least one projection indicates NOT intersecting
        logic=false;
    end
            
end

%% demo
% crv1=nrbtestcrv;
% 
% pt=rand(10,2)*5;
% [ vert,center,sa,id ] = OBB_2D( pt );
% vert(:,end+1)=vert(:,1);
% figure;
% plot(pt(:,1),pt(:,2),'b*','MarkerSize',6);
% hold on;
% plot(vert(1,:),vert(2,:),'LineWidth',2);
% quiver(center(1),center(2),sa(1,1)*5,sa(2,1)*5);
% quiver(center(1),center(2),sa(1,2)*5,sa(2,2)*5);
% axis equal;
%  
% [c1,n1]=nrbcrvbox_2D(crv1);
% nrbplot(crv1,100);
% hold on;
% num=length(c1);
% for i=1:num
%     plot(c1{i}(:,1),c1{i}(:,2));
% end
%  
% k=3;
% box1=boxconstruct_2D('OBB',vert,sa);
% box2=boxconstruct_2D('POLYGON',c1{k});
% logic=colldetect_2D(box1,box2)


