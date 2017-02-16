# Pairwise 3D Registration Evaluation

[ ![License] [license-image] ] [license]

[license-image]: https://img.shields.io/badge/license-gpl-green.svg?style=flat
[license]: https://github.com/aliosciapetrelli/Pairwise3DRegistrationEvaluation/blob/master/LICENSE

Description
-----------
C++ framework for the evaluation of Pairwise 3D registration algorithms. The framework follows the methodology proposed in:

[Petrelli A., Di Stefano L., "Pairwise registration by local orientation cues", Computer Graphics Forum, 2015.](http://onlinelibrary.wiley.com/doi/10.1111/cgf.12732/epdf)

Webpage
-----------
http://www.vision.deis.unibo.it/research/78-cvlab/108-pairwiseregistrationbenchmark

Usage
-----------
To evaluate an algorithm, the user has to implement the *PairwiseRegistration()* method of the *IPairwise3DRegistrationAlgorithm* interface. The method takes two partial views (Vi,Vj) as parameters and returns the rigid motion that aligns Vj to Vi.  

Then, the evaluation is performed by the *Pairwise3DRegistrationBenchmark* component. Once the dataset and the algorithm under evaluation have been set, the *Evaluate()* method attempts to estimate the rigid motions that align all the view pairs and, finally, computes the figure of merits.

An example can be found in [Pairwise3DRegistrationEvaluation_testmain.cpp](https://github.com/aliosciapetrelli/Pairwise3DRegistrationEvaluation/blob/master/Pairwise3DRegistrationEvaluation_testmain.cpp).

Dependencies
-----------
The framework requires [VTK](http://www.vtk.org/) and [Generalized-ICP](http://www.robots.ox.ac.uk/~avsegal/generalized_icp.html).

Note: Instead of the Generalized ICP, it is possible to apply an internal implementation of the standard ICP by uncommenting the line *#define GICP_IS_NOT_INCLUDED* in [Pairwise3DRegistrationEvaluation_FineRegistration.h](https://github.com/aliosciapetrelli/Pairwise3DRegistrationEvaluation/blob/master/Pairwise3DRegistrationEvaluation_FineRegistration.h).

The code has been tested with VTK 5.10 on Windows 7 and Microsoft Visual Studio 2010.

