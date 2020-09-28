function nurl=nrlmake(coefs, knots, intervals, order)

% Make nurl curve, surface  or  volume
% 
% Calling Sequence:
%      nurl=nrlmake(coefs, knots, intervals, order)
%      nurl=nrlmake(coefs, knots, intervals)
%      nurl=nrlmake(coefs, knots)
%  
%  Input: 
%      coefs -  Cartesian coordinates and weights
%      knots - knot vectors
%      intervals - intervals of the nurls  这个输入时到底要怎么设定？
%      order - order of the nurls basis 
%注意coefs是按行排列坐标分量，第一行是x，第二行是y，第三行是z
%  Output:
%      nurl - a nurl curve, surface or volume
% 

np = size(coefs);
dim = np(1);
m=numel(knots);
if nargin==2
    if iscell(knots)
        if m==2
            intervals={[0, 1], [0, 1]}; 
            order=[2, 2]; 
        else
            intervals={[0, 1], [0, 1], [0, 1]}; 
            order=[2, 2, 2]; 
        end
    else
        intervals=[0, 1]; 
        order=2; 
    end
elseif nargin==3
    if iscell(knots)
        if m==2
            order=[2, 2]; 
        else
            order=[2, 2, 2]; 
        end
    else
        order=2; 
    end
end

% Get numbers and reformulate coefs
if iscell(knots)
    if m==2
        m=[length(knots{1}), length(knots{2})];
        if (dim < 4)
             nrl.coefs = repmat([0.0 0.0 0.0 1.0]',[1 np(2:3)]);
             nrl.coefs(1:dim,:,:,:) = coefs;  
        else
            nrl.coefs=coefs;
        end
    else
        m=[length(knots{1}), length(knots{2}), length(knots{3})];
        if (dim < 4)
             nrl.coefs = repmat([0.0 0.0 0.0 1.0]',[1 np(2:4)]);
             nrl.coefs(1:dim,:,:,:) = coefs;  
        else
            nrl.coefs=coefs;
        end
    end
else
      if (dim < 4)
            nrl.coefs = repmat([0.0 0.0 0.0 1.0]',[1 np(2)]);
            nrl.coefs(1:dim,:) = coefs;  
      else
            nrl.coefs=coefs;
      end
end

% Make nurls curve, surface  or  volume
nurl.form='L-NURL';
nurl.dim=4;
nurl.number=m;
nurl.coefs=nrl.coefs;
nurl.knots=knots;
nurl.intervals=intervals;
nurl.order=order;










