function nrledgeplot(nurl, subd, clr)

% Plot the edges of a NURL surface or volume
% 
% Calling Sequence:
% 
%   nrledgeplot (nurl, subd)
%   nrledgeplot (nurl, subd, clr)
% 
%  Input:
% 
%    nurl - a NURL surface or volume
%    subd - number of points
%    clr - color of the edges (default value is 'k')
%

if nargin==2
    clr='k';
end

tf = ishold;
hold on;
if (iscell (nurl.knots))
      if (numel (nurl.knots) == 3)          
          srfs = nrlextract(nurl);
          for j=1:6              
              nrledgeplot(srfs(j), subd, clr);
          end
      else
          crvs = nrlextract(nurl);
          for j=1:4
              p = nrleval (crvs(j), linspace (0, 1, subd));
              plot3 (p(1,:), p(2,:), p(3,:), clr); 
          end
      end
end

if tf
    hold on;
else
    hold off;
end


%% Demo
% load torque srfs
% figure; hold on;
% for i=1:numel(srfs)
%     srfs(i).coefs(3,:)=0;
%     srfl(i)=nrb2nrl(srfs(i));
%     nrlplot(srfl(i), [10, 10],'ctrl');
% end
% axis equal;
% 
% for i=1:numel(srfl)
%     nrledgeplot(srfl(i), 10)
% end
% 
% figure; hold on;
% for i=1:numel(srfs)
%     vols(i)=nrlextrude(srfl(i), [0 0 0.1]);
%     nrlplot(vols(i), [10, 10, 10], 'ctrl');
%     nrledgeplot(vols(i), 10)
% end
% view(3); axis equal;






