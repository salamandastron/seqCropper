Prerequisites:
- RStudio installed with packages "shiny" and "shinyIncubator" installed
- Cytoscape installed with the CytoscapeRPC plugin installed
- Javascript enabled in web browser
- I am not sure but the GNU Scientific Library (GSL) may need to be installed on the computer. If this does not work, you may have to build a system-specific dynamic library file using the source codes provided in the folder /Cpp, where "commextrNP.cpp" is for the DiME algorithm itself and the "cb-signi" folder contains source for the B-score algorithm. The file to be compiled is "compare.cpp" in the subfolder "Sources" (the "-lgsl" option must be added to the linking command).

How to run:

Open the file runGUI.R in RStudio and run. Select input network file (as edgelist) and additional options, and then press "Start Algorithm".

An example network file "simNet.csv" is provided in the archive. It is an edgelist for a simulated network of 50 nodes with 3 pre-defined modules of sizes 15, 10, 5 respectively. Within-module edge probabilities are set to 0.3, 0.25 and 0.2 respectively, and inter-module and background edge probabilities are set to 0.05. "answer.txt" provides pre-defined module membership for the network.

Note that parallel computing is disabled in this version because it requires the OpenMP library which must be separately installed and configured on the computer where the program is run, and the source files will have to be re-compiled.