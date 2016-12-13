Invert shadowgraphy by Muhammad Firmansyah Kasim (University of Oxford, 2016)

This repository contains MATLAB codes to invert shadowgraphy images (or proton radiography images) to get its deflection potential of a diagnosed object.
The description about this method and some assumptions can be found on http://arxiv.org/abs/1607.04179.
It contains 3 main files:
* main_forward
* main_inverse_extended
* invert_shadowgraphy

In this code, I am using a normalised unit, where every pixel has size of 1 x 1, the distance from the object to the screen is 1, and the magnification is 1.
The normalised deflection potential (![Phi](http://sp.mfkasim.com/assets/other/invert-shadowgraphy/phi.png)) and the deflection potential in normal units 
(![Phi_N](http://sp.mfkasim.com/assets/other/invert-shadowgraphy/phiN.png)) are related by

![Phi_N = Phi * 1/L * l_pix^2 * M](http://sp.mfkasim.com/assets/other/invert-shadowgraphy/main_equation.png)

where ![L](http://sp.mfkasim.com/assets/other/invert-shadowgraphy/L.png) is the distance between the object and the screen,
![l_pix](http://sp.mfkasim.com/assets/other/invert-shadowgraphy/lpix.png) is the pixel size on the object plane, and 
![M](http://sp.mfkasim.com/assets/other/invert-shadowgraphy/M.png) is the magnification.

The codes contain other codes from:
* https://uk.mathworks.com/matlabcentral/fileexchange/56633-bounded-power-diagram
* https://uk.mathworks.com/matlabcentral/fileexchange/34428-voronoilimit
* https://uk.mathworks.com/matlabcentral/fileexchange/44385-power-diagrams
* https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html

###Getting started
######Requirements
* MATLAB R2013b or later
* [Git](https://desktop.github.com/)

######Steps
* Open cmd / bash / terminal and type the codes below

```
cd [the-directory-you-want-to-put-the-codes]
git clone https://github.com/mfkasim91/invert-shadowgraphy
cd invert-shadowgraphy
matlab
```

* Add the libraries by typing `add_libs` in MATLAB command window
* (Optional) To run the demo, type `demo` in MATLAB
* To run the inversion in your computer, you can type:

```matlab
invert_shadowgraphy('directory/to/your/proton/radiogram/image');
```

###Submitting to clusters
######Requirements
* [MCR 8.3](https://uk.mathworks.com/products/compiler/mcr.html)

######Steps
* Copy the executable in the /bin folder to your cluster
* In the script to submit the job, include this lines:

```
module load matlab/mcr/8.3
./invert_shadowgraphy /complete/path/to/your/image <options>
```

###Options
* `fdir_out`: the directory to put the results (default: pwd)
* `num_sites`: number of sites (default: min(100000, number of pixels))
* `source_map`:
  * 0 - uniform with the same size as the file (default: 0)
  * positive number - using tvdenoise of the input image with lambda = source_map argument
  * string - filename of the source
* `algorithm`: optimization algorithm, 'quasi-newton' or 'lbfgs' (default: 'lbfgs')
* `num_workers`: number of workers in the distributed computing (default: 12)
