### Renan Liupekevicius Carnielli [r.liupekevicius.carnielli@tue.nl](mailto::r.liupekevicius.carnielli@tue.nl)

Data Availability of the article "Equivalent continuum for viscoelastic metamaterials"


-*COMSOL DISPERSION AND LOCAL RESONANT MODES* is the folder with the COMSOL models to compute the local resonance modes (Figure 3); and the blue crosses DNS points of Figures 4, 5 and 6.

-*COMSOL FINITE DNS LRAM* is the folder with the COMSOL models that compute the DNS points of the transmission problem (Figure7) and the two-dimensional beam (figure 8).

-*COMSOL HOM* is the folder with the comsol models that compute the homogenised transmission problem, Figure 7 solid curves, and the two-dimensional beam, Figure 8 solid curves.



## How does 'Walrus' software works?

1- run `start.m` to include the the path to the folders *fun*, *tensorlab* and *meshes*.

2- export comsol 5.4 mesh file `.mphtxt` of your favorite geometry to the folder *meshes*, or use the mesh file provided.

3- FATAL: remember to rename the exported comsol mesh `.mphtxt` to `.txt`.

4- the files `computedworkspace_Liang.mat` and `computedworkspace_Liang_willis.mat` are the output workspaces of `main.m` for: the symmetric unit cell, figure 3(a); and the asymmetric unit cell, figure 7, respectively.

## `main.m` file description

1- `main.m` computes the homogenised material coefficients of equations (33) and equation (34).

2-  The input for `main.m` are: the mesh file in the folder *meshes* and the material properties of the saturating fluid which can be defined within the script `main.m`.

3-  The appropriate mesh is created by exporting a text file `.mphtxt` from comsol 5.4, see 'Note' 2nd and 3rd item.

## `effective_material_properties.m` file description
This script computes the frequency-domain material properties of equation (37) assuming you ran `main.m`, or loaded the provided workspaces.

## *fun* folder
Developed functions to assist the computations in `main.m`.

## *tensorlab* folder
Collection of functions and classes that define a tensor object.
