Fluorescence Imaging Toolbox for Image Systems Engineering Toolbox
(fiToolbox for ISET)

This toolbox serves as an extension of the Image Systems Engineering Toolbox 
(ISET) which provides a set of functions to simulate scenes with fluorescent 
components. The fiToolbox also implements several reflectance and
fluorescence separation and estimation algorithms. A copy of ISET is not
needed to run these algorithms.


If you use this code in your work, please cite

Henryk Blasinski, Joyce Farrell and Brian Wandell - "Simultaneous reflectance
and fluorescence spectra estimation" 2016 
(https://arxiv.org/abs/1605.04243)


0. License
----------

The code is provided as is. You are free to use and modify the code in 
non-commercial and research applications. 
If you are interested in commercial applications, please contact us as the 
method and apparatus is a subject of a US Patent 20,160,116,410.


1. Installation
---------------

To install the fiToolbox run fiToolboxInit function. This function adds all
the relevant directories to MATLAB path. Make sure that the additional
dependencies (see below) are on the MATLAB path.



2. Dependencies (required)
--------------------------

a. Image Systems Engineering Toolbox (ISET/isetcam) (https://github.com/ISET/isetcam).
b. cvx Convex Optimization toolbox for Matlab (www.cvxr.com).
c. Image Systems Engineering Toolbox requires the MATLAB Statistics and Machine Learning Toolbox
d. Image Processing Toolbox MATLAB



3. Dependencies (optional)
--------------------------

The following code is used by some data plotting scripts and are not required
to simulate or analyze reflectance and fluorescence.

a. Barwitherr function from Matlab File exchange. 
   (http://www.mathworks.com/matlabcentral/fileexchange/30639-barwitherr-errors-varargin).
b. Code accompanying the Computational Colour Science using Matlab from
   Matlab File exchange.
   (http://www.mathworks.com/matlabcentral/fileexchange/40640-computational-colour-science-using-matlab-2e)


4. Sample data
--------------

The sample inputs to our algorithms as well as the results of the analyses 
can be downloaded from Stanford's Digital Repository: 
http://purl.stanford.edu/xc528jd5098
The zip file downloaded from the above side contains two directories: 'data'
with raw input data from real captures and simulations and 'results' storing
algorithm outputs. Please place these two folders directly in the fiToolbox
root folder.




