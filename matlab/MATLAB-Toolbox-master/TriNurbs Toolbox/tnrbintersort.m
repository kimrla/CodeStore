function [sq, id]=tnrbintersort(tnrb1, ed1, pts1, p2t1)

% tnrbintersort: Sort a portion or all of the intersections of two surfaces
% 
% Calling Sequences:
% 
%       [sq, id]=tnrbintersort(tnrb1, ed1, p2t1)
% 
% INPUTS:
% 
%       tnrb1 - Triangular representation of a nurbs surface.
% 
%       ed1-  Indexes of an edge on tnrb1.
% 
%       pts1 - Approximated parametric intersections on a surface.       
% 
%       p2t1 - The relations from points to triangles.
%
% OUTPUT: 
%
%       sq  -  The sequence of index of ed1 of the intersections.
% 
%       id  -  A logical vector indictes the sorted intersections.
%

% Prepare vectors
ne=length(ed1); id=true(ne,1); sq=zeros(ne,1);

% Get two arbitrary edges and a related triangle with intersections
es=1;
tri=tnrbedge2tri(p2t1, ed1(es,:)); 
t2e=tnrbtri2edge(tnrb1, ed1, tri); 
cs=length(t2e);
nt=length(tri);
if cs==2 && nt==1
    trg=sort(RemDuplicate([ed1(t2e(1),:), ed1(t2e(2),:)]'))'; 
    tral=sort(tnrb1.delaunay(tri,:),2); 
    td=DistanceMatrix(trg, tral); 
    tri1=tri(td==0); 
    sq(1:cs)=t2e;
    id(sq(1:cs))=false;
    t=cs;
elseif cs==3 && nt==2
    te=t2e(t2e~=es);
    t2e=[te(1); es; te(2)]; 
    trg=sort(RemDuplicate([ed1(t2e(2),:), ed1(t2e(3),:)]'))'; 
    tral=sort(tnrb1.delaunay(tri,:),2); 
    td=DistanceMatrix(trg, tral); 
    tri1=tri(td==0); 
    sq(1:cs)=t2e;
    id(sq(1:cs))=false;
    t=cs;
elseif cs~=1 && nt>=1
    t2e=tnrbtri2dges(tnrb1, ed1, tri(1)); 
    t2e=sortnodes(pts1, t2e); 
    tri1=tri(1); 
    sq(1:length(t2e))=t2e; 
    id(sq(1:length(t2e)))=false; 
    t=length(t2e); 
    if nt==2
        t2e=tnrbtri2dges(tnrb1, ed1, tri(2)); 
        t2e=sortnodes(pts1, t2e); 
        if t2e(1)~=sq(t)
            t2e=t2e(length(t2e):-1:1);
            t2e=t2e(t2e~=sq(t));
            tri1=tri(2); 
            sq(t+1:t+length(t2e))=t2e; 
            id(sq(t+1:t+length(t2e)))=false; 
            t=t+length(t2e); 
        end
    end
elseif cs==1
    warning('The interection on edge [%i, %i] may be alone!', ed1(es,1), ed1(es,2));
end

% Sort the sequences of other intersection points for cs = 2
if cs>1
    while ~isempty(find(id, 1)) && t<2*ne
        % Get the next triangle
        tri=tnrbedge2tri(p2t1, ed1(sq(t),:));
        tri=tri(tri~=tri1);
        if isempty(tri)
            break;
        end
        tri1=tri(1);

        % Get the next intersected edge    
        t2e=tnrbtri2dges(tnrb1, ed1, tri(1));
        if length(t2e)>2
            t2e=sortnodes(pts1, t2e); 
            if t2e(1)~=sq(t)
                t2e=t2e(length(t2e):-1:1);
            end
        end
        t2e=t2e(t2e~=sq(t));
        if isempty(t2e)
            break;
        end
        sq(t+1:t+length(t2e))=t2e; 
        id(sq(t+1:t+length(t2e)))=false; 
        t=t+length(t2e); 
    end
    
    % Reverse the sequence
    sq(1:t)=sq(t:-1:1);
    tri=tnrbedge2tri(p2t1, ed1(sq(t),:)); 
    trg=sort(RemDuplicate([ed1(sq(t),:), ed1(sq(t-1),:)]'))'; 
    tral=sort(tnrb1.delaunay(tri,:),2); 
    td=DistanceMatrix(trg, tral); 
    tri1=tri(td==0); 

    while ~isempty(find(id, 1)) && t<2*ne
        % Get the next triangle
        tri=tnrbedge2tri(p2t1, ed1(sq(t),:));
        tri=tri(tri~=tri1);        
        if isempty(tri)
            break;
        end
        tri1=tri(1);

        % Get the next intersected edge    
        t2e=tnrbtri2dges(tnrb1, ed1, tri(1));
        if length(t2e)>2
            t2e=sortnodes(pts1, t2e); 
            if t2e(1)~=sq(t)
                t2e=t2e(length(t2e):-1:1);
            end
        end
        t2e=t2e(t2e~=sq(t));
        if isempty(t2e)
            break;
        end
        sq(t+1:t+length(t2e))=t2e; 
        id(sq(t+1:t+length(t2e)))=false; 
        t=t+length(t2e); 
    end
    sq(t+1:end)=[];
else
    sq=[];
end


%% demo
% % The mesh seed length (h0)
% h0=1.0;
% 
% % Create a nurbs sphere
% center=[9,5,5];
% circ=nrbcirc(4, center, 0, pi);
% srf1=nrbrevolve(circ, center, [1,0,0], 2*pi);
% srf2=nrbtestsrf;
% 
% % Transform a nurbs surface into triangular representation
% tnrb1=nrb2tri(srf1, h0);
% tnrb2=nrb2tri(srf2, h0);
% 
% % The nearest points of the two surfaces
% tol=max([tnrb1.seeds(1), tnrb2.seeds(1)]);
% [p1, p2, d]=nearpnts(tnrb1.points, tnrb2.points, tol);
% 
% % Get the edges of a tri-nurbs surface that intersected with another tri-nurbs surface
% p2t1=tnrbpts2tri(tnrb1);
% p2t2=tnrbpts2tri(tnrb2);
% 
% % Intersections of two tri-nurbs surfaces
% [ed1, tri2, pts1, pts2, pnts1, pnts2]=tnrbintertri(tnrb1, tnrb2, p2t1, p2t2, p1, p2);
% 
% % Sort a portion or all of the intersections of two surfaces
% [sq, id]=tnrbintersort(tnrb1, ed1, pts1, p2t1);
% 
% % Plot the results
% figure; hold on; 
% tnrbplot(tnrb1); 
% tnrbplot(tnrb2); 
% axis equal; view(3); 
% title('Geometric grid'); 
% plot3(pnts1(:,1), pnts1(:,2), pnts1(:,3), 'ro'); 
% plot3(pnts2(:,1), pnts2(:,2), pnts2(:,3), 'r*'); 
% 
% figure; hold on;
% triplot(tnrb1.delaunay, tnrb1.nodes(:,1), tnrb1.nodes(:,2));  
% % triplot(tnrb1.delaunay(tri,:), tnrb1.nodes(:,1), tnrb1.nodes(:,2), 'r');  
% plot(pts1(:,1), pts1(:,2), 'b.', 'MarkerSize', 13); 
% plot(pts1(sq,1), pts1(sq,2), 'r'); 
% % for i=1:length(t2e)
% %     ei=ed1(t2e(i),:);
% %     plot(tnrb1.nodes(ei(:),1), tnrb1.nodes(ei(:),2), 'g');  
% % end
% title('Parametric mesh of surface 1');  
% axis equal; 
% 
% figure; hold on; 
% triplot(tnrb2.delaunay, tnrb2.nodes(:,1), tnrb2.nodes(:,2)); 
% plot(pts2(:,1), pts2(:,2), 'b.', 'MarkerSize', 13); 
% title('Parametric mesh of surface 2'); 
% axis equal;








