This code is a MATLAB version which related to the paper "A New Calibration Technique for Multi-Camera Systems of Limited Overlapping Field-of-Views" that we published in IROS 2017.
 
When you clone this repository, you can simply run 'TestCheckerboardDetection.m' for quick demonstration. I strongly recommend you to run this code under GUI supported MATLAB environment. 

To use this code for calibration, you can print the pattern, take some sample images and run our code for detecting checkerboard points. For geometry estimation, I recommend you to use Zhang's method. 

If you want to replace the pattern, you have to rebuild your own 'PatternInfo.mat'.

If you use this code in your work, you have to cite our paper.

The high-performance version is implemented under C++ and it is protected by a PCT patent. For commercial usage, you need to be careful with the legal issues.
