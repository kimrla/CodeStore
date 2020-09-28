%% Global variables
function GlobalVariable(handles)
    global g_strPlane;              % Drawing view: 3D View, X-Y Plane, Y-Z Plane, Z-X Plane
    global g_bDrawControlPoint;     % Display the poles of nurbs node
    global g_IJLIncrementalMode;    % Circle programming IJK incremental mode
    global g_nScatterPrecision;     % Curve scatter precision
    global g_vecColor;              % Figure color: line point, arc point, nurbs point, control point
    
    % Initialization
    g_strPlane = handles.metricdata.Plane;
    g_bDrawControlPoint = handles.metricdata.DisplayCP;
    g_IJLIncrementalMode = handles.metricdata.IJKIncrementalMode;
    g_nScatterPrecision = handles.metricdata.ScatterPrecision;
    g_vecColor{1} = handles.metricdata.firstFileColor;
    g_vecColor{2} = handles.metricdata.secondFileColor;
end