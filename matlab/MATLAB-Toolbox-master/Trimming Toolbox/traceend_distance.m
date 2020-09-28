function lg=traceend_distance(d,tol)
% Ending condition 4: When the 2 patches are parallel, resulting to end the
% process of tracing method. (Direction Lost).

lg=1; 
if any(d>tol) % size(d)=[1,2]
    lg(d>tol)=-1;
end

end