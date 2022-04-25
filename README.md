# Nozzle-Design-MoC
This repository presents the files for the generation of the nozzle using the Method of Characteristics (MoC). MoC is a classic technique of designing the contour of the supersonic nozzle for the isentropic flow within the duct. The MoC is the direct design technique that provides a graphical solution of a flowfield, developed by Ludwig Prandtl ann Adolf Busemann. The basis of the MoC lies in attaining the contour of a nozzle by the expansion of the steady supersonic flow through Mach waves such that a supersonic, shockfree, and uniform flow is generated at the nozzle exit.

In this repository, the MoC method has been applied to attain a contour for the minimum length nozzle (MLN) in two-dimensions, programmed using MATLAB. The mathematical derivation and equations of the MoC for minimum length nozzle are available in several notable texts on fluid dynamics and compressible flows such as Modern Compressible Flow by John D. Anderson, Elements of Gas Dynamics by Liepmann and Roshko, etc.

## Inputs ##
At the beginning of the program, the user is prompted to enter the following inputs-
* Design Mach, Me
* Ratio of Specific heats, γ
* Number of Characteristics, n
* X-coordinate of sharp corner at the throat, x0
* Y-coordinate of sharp corner at the throat, y0 

The above inputs by the user are used by the program to solve the flow field using marching method. That is, the solution of the flowfield is obtained by the moving from point to point within the defined number of charcateristics, calculating the number of flow properties and direction at each node. The users who are interested in the equations employed to solve the flowfield are referred to Chapter 11 of Modern Compressible Flow by John D. Anderson. 
 
## Outputs ##

Once the program solves teh flow field, teh following output is generated-
* Flow properties (Mach number, Mach angle, deflection angle, Pradntl Meyer function) at each of the nodes.
* Coordinates of each node point in 2 dimensions and subsequently plot of the nozzle contour.
* File with data points of the MLN contour to be used as an input for Solidworks.

## Example ##
In this example, the MoC program is used to compute and graph the contour of a two-dimensional MLN for the expansion of air to a design exit Mach number, Me=2.2. The ratio of specific heat, γ =1.4 and there are 50 characteristics defined to generate the contour. The generated contour is shown in the figure below.

![Git-MoC-MLN](https://user-images.githubusercontent.com/61012294/165047236-cd0a830e-b932-4071-940b-ef37546b49ac.png)
