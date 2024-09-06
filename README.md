# MAGBOT

## Experimental Setup: 
The experimental system requires an LED platform to be laid flat on the ground, with a camera placed horizontally above the platform. The magbot moves on the LED platform, and the camera captures its positional information via a light source on top of the magbot. MATLAB is used to control the magbot swarm to perform specific experimental tasks through the following code.
After the experimental setup is completed, MATLAB's built-in Camera Calibration Toolbox is used to calibrate the camera, and the resulting parameters are stored in the file cameraParams_MAGBOT.mat. Additionally, the WindowAPI file is used to project output spots onto the LED platform, and this file can be downloaded from the following URL.https://ww2.mathworks.cn/matlabcentral/fileexchange/31437-windowapi
![figure1](https://github.com/user-attachments/assets/93126bfa-ace0-4d25-899a-2243065d24bf)


## Auto_Coordinated_internal_motion
This code is used in the Coordinated internal motion experiment to automate the control of two smaller green clusters, rotating counterclockwise around the larger blue central cluster for half a circle.

## Auto_Infiltration
This code is used in the Infiltration experiment to automate the control of small clusters moving through narrow channels to the target area.

## HCI_Smart_morphing
This code is used in the Smart Morphing experiment to manually control small clusters to complete the task objectives.

## HCI_Infiltration
This code is used in the Infiltration experiment to manually control small clusters to complete the task objectives.

## Simulation
This code is used to perform the dynamic simulation of the Magbot system. See Section 4 in the Supplementary Information for the definition of parameters.  
