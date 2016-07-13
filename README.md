Invert shadowgraphy by Muhammad Firmansyah Kasim (University of Oxford, 2016)

This repository contains MATLAB codes to invert shadowgraphy images (or proton radiography images) to get its deflection potential of a diagnosed object.
The description about this method and some assumptions can be found on ??? (TBC).
It contains two main files:
* main_forward
* main_inverse_extended (coming soon)
* invert_shadowgraphy (coming soon)

In this code, I am using a normalised unit, where every pixel has size of 1 x 1, the distance from the object to the screen is 1, and the magnification is 1.
The normalised deflection potential ($\Phi$) and the deflection potential in normal units ($Phi_N$) are related by
$$
\Phi_N = \Phi * 1/L * l_{pix}^2 * M,
$$
where $L$ is the distance between the object and the screen, $l_{pix}$ is the pixel size on the object plane, and $M$ is the magnification.
