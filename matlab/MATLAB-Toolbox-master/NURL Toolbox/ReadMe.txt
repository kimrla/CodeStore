This file presents brief introductions of a collection of routines for the creation and operation of Non-Uniform Rational Lagrange (NURL) geometries. The package includes 3 sub-packages: (1) a basic NURL package. See A1 to A8. Due to the equivalence between NURL and NURBS geometries, and for convenience of those who are familiar with the NURBS Toolbox, the data structure and the inputs and outputs of this package are coded to be similar as the NURBS package. Some NURBS routines can be or have been directly used in the NURL package. Therefore, one need to set the NURBS package in roots before using the NURL packages. (2) A NURL curve package. See A9. This package includes routines for creation and operation of plane NURL curves such as curve creation, cutting, etc. (3) A graphic user interface (GUI) based on the NURL curve package. See A10. This package can be used to create plane Coons surfaces using plane curves. Most routines include detailed descriptions and a few demos. The way of using the routines should be understood through the following description and reading the descriptions included in the codes and the demos therein. The NURL packages can be downloaded from this link: https://sourceforge.net/projects/nurl/.

A1. Basic operations for NURL curves, surfaces and volumes

nrlmake: Construct the NURL structure through given the control points and the knots.
nrlkntins: Insert a single or multiple knots into a NURL curve, surface or volume.
nrlintins: Insert a single or multiple intervals into a NURL curve, surface or volume.
nrldegelev: Elevate the degree of the NURL curve, surface or volume.
nrleval: Evaluate a NURL curve, surface or volume at parametric points.
nrldeval: Evaluation of the derivative and second derivatives of NURL curve, surface or volume.
nrldseval: Evaluation of high order derivatives of NURL curve, surface or volume.
nrlgeval: Evaluate a NURL at parametric points by interpolation matrix.
nrlgdeval: Evaluation of first and second derivatives of NURL by interpolation matrix.

A2. Operations for constructing NURL curves and surfaces

nrltform: Apply transformation matrix to the NURL.
nrlreverse: Reverse the evaluation directions of a NURL geometry.
nrltransp: Transpose a NURL surface, by swapping U and V directions.
nrlpermute: Rearrange the directions of a NURL volume or surface.
nrlmeasure: Measure the length of a NURL curve.
nrlline: Construct a straight line.
nrlcirc: Construct a circular arc.
nrlrect: Construct NURL representation of a rectangular curve.
nrl4surf: Constructs a NURL bilinear surface.
nrlcylind: Construct a cylinder or cylindrical patch.
nrlellip: Create a plane ellipse nurl curve.
nrlconic: Create a nurl conic arc (less than pi).
spiralcylind: Create a interpolated cylindrical spiral curve.
nrlsphere: Create a nurl sphere.
nrlextract: construct NURL curves by extracting the boundaries of a NURL surface, or NURL surfaces by extracting the boundary of a NURL volume.
nrlextrude: Construct a NURL surface by extruding a NURL curve, or construct a NURL volume by extruding a NURL surface.
nrlrevolve: Construct a NURL surface by revolving a NURL curve, or construct a NURL volume by revolving a NURL surface.
nrlruled: Construct a ruled surface between two NURL curves.
nrlcoons: Construction of a Coons patch.
crvgintvarcder: Evaluation of the first and second derivatives of NURL curve in arc length coordinates.
curvature: Get the curvature of a curve.
nrlcvttors: Get the curvature and torsion of a nurl curve.
nrlglue: Glue two nurl curves together.
nrlaxis: Get the axises of a nurl geometry.
nrlxyz: Get the x, y, z coordinates of a nurl surface.
nrltestcrv: Constructs a simple test curve.
nrltestsrf: Constructs a simple test surface.
torque: Test NURBS surfaces.

A3. Plot and import

nrlplot: Plot a NURL curve or surface, or the boundary of a NURL volume.
nrblplot: Plot a NURL curve or surface, or the boundary of a NURL volume in the same manner as nrbplot of NURBS toolbox.
nrlctrlplot: Plot a NURL entity along with its control points.
nrlcrvplot: Plot nurl curves and return their figure handles.
nrlsrfplot: Plot a nurl surface with highlighted edges and return its figure handle.
nrledgeplot: Plot the edges of a NURL surface or volume.
nrlaxisplot: Plot the axises of a nurl geometry.
nrb2nrl: Transform nurbs curve, surface or volume into nurl.

A4. The NURL for analysis

nrlgintveval: Evaluation weights at parametric points. 
nrlgintvdeval: Evaluation matrices of first and second derivatives of NURL curve, surface or volume.
srfgintvplaneder: Evaluation weights of the first derivatives of NURL plane surface in Cartesion coordinates.
VibBarL: Vibration analysis of bars by NURL.
VibBarH: Vibration analysis of bars by NURH.

A5. NURL functions for triangles

nrl2trg: Transform a nurl triangular patch for subsequent manipulations.
nrltrgcoons: Construction of a triangular Coons patch.
srftrgplaneder: Evaluation weights of the first derivatives of a NURL triangular plane surface in Cartesian coordinates.
nrltrgdeval: Evaluate the derivatives of an triangle with respect to area coordinates.
nrltrgmat: Get first derivatives of the nurl basis matrix in area coordinates of a triangle.
trgsrfdirect: Rearrange the direction (axis) of a triangular nurl patch for subsequent manipulations. 

A6. Non-Uniform Lagrange (NUL) functions

spanweak: Get the spans and indexes of evaluation points like the findspan function of NURBS.
nurlbasis: Get non-uniform Langrane basis.
nurhbasis: Get the basis of non-uniform Hermite basis.
DeriveBasis: Deriving the nurl basis by symbolic computation.
nulpts2crv: Get a nul curve through arbitrary points.
nulkntins: Insert knots into a NUL.
nulintins: Insert intervals into a NUL.
nuldegelev: Degree elevation of NUL.
nulintvdeval: Evaluate nul curves or their derivatives at an interval.
nuldeval: Evaluate nul curve, surface, volume or their derivatives.
nulintvmat: Weighting coefficient matrix of Galerkin interpolation or its derivatives at an interval for NUL basis.
nuhintvmat: Weighting coefficient matrix of Galerkin interpolation or its derivatives at an interval for NUH basis.

A7. Related routines for geometries and analysis

DistanceMatrix: Forms the distance matrix of two sets of points in Rs.
FindUnique: Find unique elements for a (n*d) vector matrix.
RemDuplicate: Remove duplicate elements for a (n*d) vector matrix.
GaussLobattoR: Solve roots and weights of Gauss-Lobatto quadrature.
GaussLobattoQ: Roots and weights of Gauss-Lobatto quadrature with intervales.
GaussR: Solve the nodes and weights of Gauss integration method.
LobattoChebyshev: Generate Lobatto Chebyshev nodes.
planeline: Create plane lines and get relations of plane lines through analytical expressions.
planepls: Get relations between planes and plane lines.
countknts: Count the number of knots in each intervals.
checktt: Check whether a cell array is with row vectors.
SolveEig: Solve eigenvalues.
Weighting: Get weighting coefficient matrix of DQM.

A8. Vector and Transformation Utilities

The NURL routines for this part are the same as the NURBS, so please refer to the NURBS toolbox. 

A9. NURL for plane curves

aintsctpnts: Get approximated points of intersection of two curves using midpoint finite difference method.
bertrand: Get double curves of a curve with a distance (known as Bertrand curves).
BisectVector: Get the bisector vector of two plane vectors.
cmpoints: Get index of the common points of two curves.
crvnearpnt: Get the nearest distance of a point to a curve.
crvrotate: Apply a rotation to a plane curve around a point.
intsctpnt: Solve a intersection point of two curves using Newton-Raphson's method.
intsctpnts: Solve points of intersection of two curves.
isnrlline: Check whether a nurl curve is a straight line.
mfdm: Get derivatives of a function using midpoint finite difference method.
neartangent: Get tangent vector of an end of a line that is nearest to a point
nrl2nrb: Transform a nurl curve to a nurbs curve.
nrl3pntsarc: Create a nurl arc (<= pi) by three points in sequence.
nrl3ptscirc: Create a nurl circle by three points.
nrlanalytical: Create a curve using analytical functions.
nrlcircenter: Get the centre of a nurl circle.
nrlcrvextract: Extract the two end points of a curve and get the derivatives of the curve at the two points.
nrlcuts: Cut curves by one curve.
nrlfillet: Filleting two straight lines.
nrlglues: Glue several curves to form a single curve.
nrlpolygon: Create a nurl polygon in x-y plane.
nrlquadr: Create a nurl quadrangle formed by 4 points.
nrlspline: Create a nurl curve using control points of B-splines.
nrlsplit: Split a nurl curve into two separate curves.
nrlsplits: Split a nurl curve into several separate curves.
nrltriangle: Create a nurl triangle formed by three points.
plangle: Get the angle of a plane vector with respect to x-coordinate.
plcrvmirror: Get the mirror of plane curve(s).
plcrvtranslat: Apply a translation to a plane curve.
remplines: Remove lines that are actually points.
sort4nodes: Sort 4 nodes for the construction of a NURL bilinear surface.
tangentarc: Get the tangent arcs of two straight lines.
tangentarcs: Get the tangent arcs of two straight lines (has more outputs compared with tangentarc).

A10. GUI for plane geometries

PlaneGeom: GUI for plane geometries.
edge2coons: Get edges for making coons surfaces from arbitrary edges.
nrlsrfsdirection: Transform the directions of plane surfaces to have a common normal vector.
CreatCoons: Create NURL Coons surfaces from NURL or NURBS curves.
brim: A test NURL geometry for CreatCoons.
crown: A test NURL geometry for CreatCoons.
torquearm: A test NURBS geometry for CreatCoons.
