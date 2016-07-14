Invert shadowgraphy by Muhammad Firmansyah Kasim (University of Oxford, 2016)

This repository contains MATLAB codes to invert shadowgraphy images (or proton radiography images) to get its deflection potential of a diagnosed object.
The description about this method and some assumptions can be found on ??? (TBC).
It contains two main files:
* main_forward
* main_inverse_extended (coming soon)
* invert_shadowgraphy (coming soon)

In this code, I am using a normalised unit, where every pixel has size of 1 x 1, the distance from the object to the screen is 1, and the magnification is 1.
The normalised deflection potential (Phi) and the deflection potential in normal units (Phi_N) are related by

![(Phi_N = Phi * 1/L * lpix^2 * M)](http://www.sciweavers.org/tex2img.php?eq=%5CPhi_N%20%3D%20%5CPhi%20%2A%201%2FL%20%2A%20l_%7Bpix%7D%5E2%20%2A%20M%2C&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

where L is the distance between the object and the screen, lpix is the pixel size on the object plane, and M is the magnification.

The codes contain other codes from:
* https://uk.mathworks.com/matlabcentral/fileexchange/56633-bounded-power-diagram
* https://uk.mathworks.com/matlabcentral/fileexchange/34428-voronoilimit
* https://uk.mathworks.com/matlabcentral/fileexchange/44385-power-diagrams
* https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html

License for non-commercial and non-military uses.
