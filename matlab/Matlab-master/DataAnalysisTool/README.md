# Data Analysis Tool Introduction

## Tool purpose and generation

This tool is used to analyze the collected track data.

Use the command `guide` to open `DataAnalysisTool.fig`, you can view and modify the interface layout and function buttons of the tool.

Use the command `mcc -m DataAnalysisTool.m` to compile and generate the executable file `DataAnalysisTool.exe`.

## Parameters

### Sampling period

The time interval at which the track data was collected.

### Data format

The storage format of the collected data in each line of the file, where `%f` represents a numerical value. If the data contains specific characters, such as `x`, `y`, etc., you can add `%c` to the data format. If the data contains punctuation marks, such as `,`, `;`, etc., you can directly add the corresponding punctuation marks to the data format. Then you can guarantee that the data points stored in the file are read correctly.

### Valid data columns

The `%f` in the data format represents a numerical value. You can specify whether a column corresponding to a certain `%f` is valid. If it is invalid, the column value will be discarded when filtering the data.

### Delete redundant data at the beginning and end

This option can remove redundant data if the collected data contains the same row data at the beginning and end.

## Drawing options

### Single axis mode

In this mode, you can only draw the time-speed, time-acceleration, and time-jerk of a single axis. You cannot draw the track and distance-speed.

### Draw track

For three dimensional or less track data, you can draw track. Track data exceeding three dimensions does not support drawing track.

The rest of the drawing options can be identified based on the parameter names, so I won't repeat them here.

## Cut options

A piece of collected data may be very large, and sometimes it is necessary to intercept a part of the data for analysis. In this case, you can use the cut option to achieve the purpose. According to the time-speed figure, determine the start time and end time to be cut. Only three sets of data can be cut at the same time.

## Importing Files

Select the files to be analyzed. After importing, you can view the imported files in the data files list in the lower left corner.

## Drawing Comparison

Import the files to be analyzed, set the corresponding parameters, and click this button to draw the corresponding figures.

## Perform Cut

Determine the cut ranges, enter them into the cut options, and click this button to draw the corresponding figures of the cut parts.

## Clear Figure

Close all figures previously drawn by `Drawing Comparison` and `Perform Cut`.

## Close Tool

`Clear Figure`, and then close the tool.