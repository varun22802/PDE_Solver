# MA5720_PDE_Solver
The function takes inputs of the Poisson Problem PDE and outputs the value of the function u inside the mentioned domain.

When main.cpp is executed, it asks for multiple inputs. It assumes A,B,C,D,E as the coefficients to d2u/dx2, d2u/dy2, du/dx, du/dy, constant term respectively. ax,bx,ay,by correspond to the diagonally opposite points that define the rectangular domain. nx,ny is the number of partitions into which x and y axis of the domain are discretized respectively. It takes the inputs f1,f2,g1,g2 as the 4 boundary conditions of the pde. The output is the value of the function u at all the points inside the mentioned domain with the given discretization.
