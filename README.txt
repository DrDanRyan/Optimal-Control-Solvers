This project is a set of optimal control solver routines for Matlab.  
Please see the OCSolverManual.pdf file for a complete description.

This project relies on the ode solver odevr7 created by Larry Shampine and 
freely available at http://faculty.smu.edu/shampine/current.html. This solver
is significantly faster at solving odes to very tight tolerances than any of 
the basic Matlab routines.

Also, the griddedInterpolant class is used extensively, which requires Matlab 2011b
or later.