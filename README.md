# Heat-diffusion-fins-radiator

## Introduction
In this project we are going to study the diffusion of heat on a fins' radiator. The problem is formulated with the relevant math equations in the document **simulation.pdf** as it is more convenient to have a separate file for LaTeX equations. In this **Readme.md** we are going to pactically describe how run the simulation and get those results.
![Demo](https://user-images.githubusercontent.com/16581022/34326109-4645c7ba-e89d-11e7-9b0f-33a6615bc7df.gif)
![Demo](https://user-images.githubusercontent.com/16581022/34326126-9bcd3e16-e89d-11e7-8148-3dfce9371c9e.gif)

## Required Softwares
* **C++** : I personally use  Xcode on macOS 10.12, but you can download Code::Blocks as an IDE and a compiler like GCC. 
http://www.codeblocks.org/downloads
* **Freefem++** is a Mesh generating software : you can have all the relevant information on the following website : http://www.freefem.org/
* **Gnuplot** is a command line graphical program to generate two and three dimensional plots of functions data and data fits.
https://sourceforge.net/projects/gnuplot/

## Folders and File Description

Before going to much in detail, if you want to run the simulation, you want all the files of this repository to be in one same working directory. I have separated them in folders to be more concise.

* **include** contains the files : *RNM.hpp*, *RNM_op.hpp*, *RNM_opc.hpp*, *RNM_tpl.hpp*. They are required in the main program in part to interpret the output of Freefem++. They are mainly used to declare multidimensional arrays such as matrices and tensors. You don't need to look at them to successfully run the simulation. To get an idea of what these files are doing you can look at the file *Test/Ex_Utilisation_RNM.cpp*. You can find additional information at the following address : https://searchcode.com/file/25416135/freefem++-3.19-1/src/femlib/RNM_op.hpp

* **Mesh_Generation** : This folder have two files in it *domain.edp* that is written using the Freefem++ syntax. When compiling you create the file *D.msh* and output the following image.

![Alt text](Mesh_Generation/maillage.png?raw=false "Title")

* *sfem.hpp* produces a set of useful classes for our simulation, in particular, it produces the following:
  * R and R2 are classes to implement real numbers and two dimensional real vectors respectively.
  * Label is used to label the vertices of a triangle in the mesh. It is useful to model specific boundary conditions to know if a vertice lie in the open domain or on the boundary.
  * Vertex inherit from the classe R2 and label. Each vertice is now considered as a vector and a label, to inform about it's position on the domain.
  * Triangle is an array of pointers on three vertices and a label, it also computes some useful properties such as area and other functions used to compute the terms of the matrix in the discretize linear system. 
  * Boundary edge inform the .main.cpp program if an edge is within the domain, or on a specific boundary. This information is important to localize Dirichlet boundaries.
  * Mesh is used to model the whole Mesh of a domain by informing the .main.cpp program about
    * the number of vertices (nv), number of triangles (nt), number of boundary edges (neb).
    * the array of vertices.
    * the array of triangles.
    * the array of boundary edges.
* *GC.hpp* implement the conjugate gradient method to solve a well conditionned linear system. An example of how this method works is given in the file *Test/Ex_Utilisation_GC.cpp*
* *main.cpp* contains all the functions and classes to execute the simulation. In particular, those are worth noticing:
  * VirtualMatrice and addMatMul modeled the matrix described in the equation (1.6) of the document **simulation.pdf**
  and compute the matrix-vector multiplication needed to solve the equation (1.5) with the conjugate gradient method
  * L2 computes the l2 norm between two vectors, it is used in the main function to evaluate the termination criterion (convergence of the temperature solution)
  * The main function initialize the temperature at t=0 and computes all the solution of the temperature until convergence.
