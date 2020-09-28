function lg=detercond(mode,varargin)
% Different determine conditions, including 'tolerance', 'iterative',
% 'mountain', 'derivative', 'cos'.

% Input:
%   mode: 'tolerance', 'iterative', 'mountain', 'derivative', 'cos'.
%   varargin: Corresponding input variances of different mode.

% Output:
%   lg: Logical scalar where 1: the iterative condition has arrived the 
%       endding; 0: Not endding which needs to be continuted to iterative. 

lg=false;

switch mode
    % The norm between the (k+1)th value and the kth value.
    case 'tolerance' % xk+1, xk=varargin{1}, varargin{2}. tol=varargin{3}
        if length(varargin)~=3
            error('The input is error!');
        end
        if norm(varargin{1}-varargin{2})<=varargin{3}
            lg=true;
        end           
    % Iterative times     
    case 'iterative' % n=varargin{1}: iterative times of the current 
        % circulation;  N=varargin{2}: default iterative times.
        if length(varargin)~=2
            error('The input is error!');
        end
        if varargin{1}>=varargin{2}
            lg=true;
        end        
    % Mountain factor    
    case 'mountain' % x=varargin{1}: factor of the current step. 
        % X=varargin{2}: factor tolerance, tol=varargin{3}
        if length(varargin)~=3
            error('The input is error!');
        end
        if abs(varargin{1}-varargin{2})<=tol
            lg=true;
        end       
    % Derivative (vector) equals 0    
    case 'derivative'
        if length(varargin)~=1
            error('The input is error!');
        end
        if norm(varargin{1})<=eps
            lg=true;
        end        
    % Cosine of the 2 vectors equals 0    
    case 'cos' % v1=varargin{1}, v2=varargin{2}. tol=varargin{3}
        if length(varargin)~=3
            error('The input is error!');
        end
        if abs(dot(varargin{1},varargin{2}))<=tol
            lg=true;
        end        
    otherwise
        error('The iterative mode is error!');
end

end