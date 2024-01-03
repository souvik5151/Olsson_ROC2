# Olsson.wl

A *Mathematica* package for deriving linear transformations of multivariate hypergeometric functions based on the method of Olsson [1]. It is based on the work :

[Olsson.wl : a Mathematica package for the computation of linear transformations of multivariable hypergeometric functions](https://arxiv.org/abs/2201.01189)

# ROC2.wl

The companion package *ROC2.wl* can find the region of convergence of double variable hypergeometric functions.


# Installation

To install the package one can copy and paste the package at the desired **location**. The path for the directory can be set as follows

    SetDirectory[Path of the location]

After setting the path, The **Olsson.wl** package can be called using the following command

    <<Olsson.wl

The comparinon package **ROC2.wl** is automatically called inside the **Olsson.wl** package.

The available commands can be found using the following command :

    ?Olsson`* 

The commands and usage of the packages are discussed in the article in great detail.


# Supplementary files 

F4_derivation.nb : This file contains derivation of analytical continuations of Appell F4 using the package. 


KdF_derivation.nb : This file contains the derivation of analytic continuations relevant to the study of Feynman integral considered in the paper.


Ac_all_combined.nb : This file contains all the combined analytic continuations relevant to the calculations of Feynman integral considered in the paper. 

transformations_HMS.nb : This file contains derivation of various transformation formulas from the book [2].



[1] [Integration of the Partial Differential Equations for the Hypergeometric Functions F1 and FD of Two and More Variables](https://pubs.aip.org/aip/jmp/article/5/3/420/230849/Integration-of-the-Partial-Differential-Equations)


[2] Srivastava, H. M., & Karlsson, P. W. (1985). Multiple Gaussian hypergeometric series.

