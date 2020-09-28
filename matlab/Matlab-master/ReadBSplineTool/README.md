# Read B-spline Tool Introduction

## Tool purpose and generation

This tool is used to read G code in the file and draw the corresponding track. Among them, the G code identification of nurbs node is G06.2, and its form is as follows

```
G06.2 P3 K0.0 X1.0 Y0.0 Z0.0 R1.0   % P is the degree of nurbs node
K0.0 X1.0 Y1.0 Z0.0 R0.5            % K is the knot vector of nurbs node
K0.0 X-1.0 Y1.0 Z0.0 R0.5           % X, Y, Z are the poles of nurbs node
K0.0 X-1.0 Y0.0 Z0.0 R1.0           % R is the weight of pole of nurbs node
K1.0
K1.0
K1.0
K1.0
```

Use the command `guide` to open `ReadBSplineTool.fig`, you can view and modify the interface layout and function buttons of the tool.

Use the command `mcc -m ReadBSplineTool.m` to compile and generate the executable file `ReadBSplineTool.exe`.

## Parameters

### View

Views when drawing track. The available views are: X-Y Plane, Y-Z Plane, Z-X Plane, and 3D View.

### Display poles of nurbs node

Whether to draw the poles of nurbs node when drawing track. If checked, it is drawn, if not checked, it is not drawn.

### Circle programming IJK incremental mode

When the G code represents a circular arc, whether the center of the circle indicated by IJK is an incremental mode or not. If it is not checked, it is an absolute value. The value of IJK directly represents the center of the circle, otherwise you need to add the start point to use as the coordinates of the center of the circle.

### Curve scatter precision

Use scattered points to replace the curve (arc or nurbs node) when drawing track. The precision here means the precision of approximating the curve with scattered points.

## Processing files

You can manually enter the file path to be analyzed, or import the file path through the button `Importing Files`.

## Drawing color

Here you can set the colors corresponding to different line shapes in the file.

## Importing Files

Select the files to be analyzed. After importing, you can view or modify the imported files in the processing files box.

## Drawing Comparison

Import the files to be analyzed, set the corresponding parameters, and then click this button to draw the corresponding track of the files.

## Clear Figures

Clear all figures previously executed by `Drawing Comparison`.

## Close Tool

Clear all figures, and then close this tool.

