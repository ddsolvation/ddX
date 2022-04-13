# The dd-story
*“A quantum chemist, an applied mathematician and a molecular dynamics expert enter in a meeting room…”*



### A long story short
#### ddCOSMO

The first aspect that should be apparent in the ddCOSMO story is that this is a multidisciplinary story, involving three different communities: applied mathematicians, quantum chemists and the classical molecular dynamics community. Therefore, the answer to the questions “what was the most rewarding part of the project” and “what was the most difficult part of the project” is unique: working together. The incredible language barrier between mathematicians and chemist has made the communication quite difficult until a common ground was found, and, after quite a few very funny reunions which ended with people looking at each other with a very puzzled face, we finally managed to fully understand each other. This certainly required a big effort from every side, but such an effort was amply rewarded not only by the results, which couldn’t have been obtained without such a collaboration, but also, or even mainly, for the fun we had working together and for the pleasure of collaborating with scientist that do not share the same education and view on problems, but can provide, every day and in every occasion, an exciting and new perspective on old and new problems.

The project has its roots in an awarded France-Berkeley-Fund proposal in 2011 which allowed the collaboration between the applied mathematicians Eric Cancès, Yvon Maday (Paris region) and Benjamin Stamm (UC Berkeley at that time). The initial plan was to use reduced order techniques in solvation models, but, as so often in life, we deviated from the original route and the ddCOSMO was developed. 

The chemistry side of this story starts in June 2012. Realizing the potentialities of the algorithm for chemical applications, Eric Cancès contacted Benedetta Mennucci to inquire whether she was interested in the development of the method, and whether she had a student who could work on the implementation and on the interface with the various quantum chemical packages. On Thursday June the 7th, Filippo Lipparini, at the time a Ph.D. student in Prof. Mennucci’s group, attended the first meeting in Paris. It was love at first sight!

A first, preliminary, pedestrian implementation was obtained in the following summer. The iterative solver was a simple Jacobi Solver, Spherical Harmonics – needed for the discretization – were generated using library calls for the trigonometric functions, no vectorization attempt of the code had been made and the stand-alone, Fortran 77 initial code was quite a mess: nevertheless, it was already at least two order of magnitude faster than the traditional COSMO implementations.

A first paper describing the algorithm was submitted to the Journal of Chemical Physics in March 2013. At the same time, the work on a more efficient implementation, and on the implementation of analytical derivatives of the ddCOSMO energy with respect to the nuclei’s positions had started. The results of this work were submitted to the Journal of Chemical Theory and Computation in April 2013. 

With the preliminary development of the ddCOSMO machinery completed, the interfacing work started. From the quantum chemistry side, a first implementation was achieved in Gaussian together with Benedetta Mennucci, first for semiempirical methods, then, in collaboration with Giovanni Scalmani, for general Quantum Mechanical models. From the classical side, ddCOSMO was implemented in Tinker for classical force field simulations first and for polarizable force fields later together with Jean-Philip Piquemal.

The results were very encouraging. We compared ddCOSMO to the fastest implementation we knew of – based on the Fast Multipole Method – and obtained a gain of two to three orders of magnitude. The comparison with traditional, quadratic scaling implementations was even more impressive, with elapsed time of several hours being reduced to a few seconds.
ddCOSMO proved to be so fast that we considered using it to provide a boundary condition in polarizable molecular dynamics simulation: taking fully into account the mutual polarization between the solvent and the solute, we could carry out MD simulations of medium-sized molecules on a standard computer node.

**A funny side note:** The first paper describing the algorithm that was submitted to the Journal of Chemical Physics in March 2013 by the three Mathematicians Cancès, Maday and Stamm received a funny review with statements such as:
*“..., however the manuscript suffers greatly from the reliance on the authors trying to write and act like mathematicians”*
*“There's a useful numerical method here, at the heart of it all, but the presentation is a convoluted mess …, if it is to be useful to (and understood by) the chemical physics community.”*
*“It's just a 3-d differential equation that they are mapping onto a boundary-element method and solving via domain decomposition. It's not rocket science, it's just linear equations -- so why present the material in a complicated (and tediously long) way?”*

#### ddPCM

#### ddLPB
