# FoldingscreenAnalysis

Analysis of TU Munich, Dietzlab - Standard Initial Folding Screen for DNA Origami
This MATLAB Tool helps determine the best folding conditons for the given design, 
    as well as folding yield, folding quality and migration distance of at the given conditions.

The Tool is composed of 3 seperate steps that can either be executed at once, or individually. 

# Usage
 Requires the following files to be present in a folder:
* gel_info.txt
* gel_image.tif

To start the full analysis:
``` IFS_analysis ```

If you need to calculate the individual lane profiles:
``` main_part1 ```
To classify the different lanes and bands:
``` main_part2 ```
For gernerating optimal folding conditions and analysis overview:
``` main_part3 ```
each substep depends on the execution of its predecessors.

# Requirements:
## External
* None

## MATLAB Toolboxes
* Image Processing Toolbox
* Curve Fitting Toolbox