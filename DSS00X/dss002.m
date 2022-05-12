%  File: dss002.m
%
   function [ux]=dss002(xl,xu,n,u)
%
%  Function dss002 computes the first derivative, u , of a
%                                                  x
%  variable u over the spatial domain xl le x le xu.
%
%  Argument list
%
%     xl      Lower boundary value of x (input)
%
%     xu      Upper boundary value of x (input)
%
%     n       Number of grid points in the x domain including the
%             boundary points (input)
%
%     u       One-dimensional array containing the values of u at
%             the n grid point points for which the derivative is
%             to be computed (input)
%
%     ux      One-dimensional array containing the numerical
%             values of the derivatives of u at the n grid points
%             (output)
%
%  Function dss002 computes the first derivative, u , of a
%                                                  x
%  variable u over the spatial domain xl le x le xu from the
%  classical three-point, second-order finite difference approxi-
%  tions
%
%                                       2
%  u1  = (1/2dx)(-3u1 + 4u2 - u3) + O(dx ) (left boundary,     (1)
%    x                                         x = xl)
%
%                                   2
%  ui  = (1/2dx)(ui+1 - ui-1) + O(dx ) (interior point,        (2)
%    x                                   x ne xl, xu)
%
%                                          2
%  un  = (1/2dx)(3un - 4un-1 + un-2) + O(dx ) (right boundary, (3)
%    x                                            x = xu)
%
%  Equations (1) to (3) apply over a grid in x with corresponding
%  values of the function u(x) represented as
%
%   u1      u2       u3         ui        un-2      un-1    un
%
%  x=xl  x=xl+dx  x=xl+2dx ... X=xi ... X=xu-2dx  x=xu-dx  x=xu
%
%  The origin of equations (1) to (3) is outlined below.
%
%  Consider the following polynomial in x of arbitrary order
%
%                                     2             3
%  u(x) = a0 + a1(x - x0) + a2(x - x0)  + a3(x - x0)  + ....   (4)
%
%  We seek the values of the coefficients a0, a1, a2, ... for a
%  particular function u(x).  If x = x0 is substituted in equation
%  (4), we have immediately a0 = u(x0).  Next, if equation (4) is
%  differentiated with respect to x,
%
%                                                   2
%  du(x)/dx = u (x) = a1 + 2a2(x - x0) + 3a3(x - x0)  + ...    (5)
%              x
%
%  Again, with x = x0, a1 = du(x0)/dx = u (x0).  Differentiation
%                                        x
%  of equation (5) in turn gives
%
%  d2u(x)/dx2 = u  (x) = 2a2 + 6a3(x - x0) + ...
%                2x
%
%  And for x = x0, a2 = u  (x0)/2f (2f = 1*2, i.e., 2 factorial).
%                        2x
%
%  We can continue this process of differentiation followed by the
%  substitution x = x0 to obtain the successive coefficients in
%  equation (4), a3, a4, ...  Finally, substitution of these co-
%  efficients in equation (4) gives
%
%                                                 2
%  u(x) = u(x0) + u (x0)(x - x0) + u  (x0)(x - x0)  +
%                  x       1f       2x       2f
%                                                              (6)
%                                3                  4
%                 U  (x0)(x - x0)  + u  (x0)(x - x0)  + ...
%                  3x       3f        4x       4f
%
%  The correspondence between equation (6) and the well-known
%  Taylor series should be clear.  Thus the expansion of a
%  function, u(x), around a neighboring point x0 in terms of u(x0)
%  and the derivatives of u(x) at x = x0 is equivalent to approxi-
%  mating u(x) near x0 by a polynomial.
%
%  Equation (6) is the starting point for the derivation of the
%  classical finite difference approximations of derivatives such
%  as the three-point formulas of equations (1), (2) and (3).  We
%  will now consider the derivation of these three-point formulas
%  in a standard format that can then be extended to higher
%  multi-point formulas in other subroutines, e.g., five-point
%  formulas in routine dss004.
%
%  Three-point formulas
%
%     (1)  Left end, point i = 1
%
%  If equation (6) is written around the points x = xl for x = xl +
%  dx and x = xl + 2dx, for which the corresponding values of u(x)
%  are u1, u2 and u3 (u1 and u2 are separated with respect to x by
%  distance dx as are u2 and u3, i.e., we assume a uniform grid
%  spacing, dx, for independent variable x)
%
%                                2            3
%  u2 = u1 + u1 ( dx) + u1  ( dx)  + u1  ( dx)  + ...          (7)
%              x  1f      2x  2f       3x  3f
%
%                                2            3
%  u3 = u1 + u1 (2dx) + u1  (2dx)  + u1  (2dx)  + ...          (8)
%              x  1f      2x  2f       3x  3f
%
%  We can now take a linear combination of equations (7) and (8)
%  by first multiplying equation (7) by a constant, a, and equa-
%  tion (8) by constant b
%
%                                  2           3
%  a(u2 = u1 + u1 ( dx) + u1  ( dx) + u1  ( dx) + ...)         (9)
%                x  1f      2x  2f      3x  3f
%
%                                  2           3
%  b(u3 = u1 + u1 (2dx) + u1  (2dx) + u1  (2dx) + ...)        (10)
%                x  1f      2x  2f      3x  3f
%
%  Constants a and b are then selected so that the coefficients of
%  the u1  terms sum to one (since we are interested in obtaining
%        x
%  a finite difference approximation for this first derivative).
%  Also, we select a and b so that the coefficients of the u1
%                                                            2x
%  terms sum to zero in order to drop out the contribution of this
%  second derivative (the basic idea is to drop out as many of the
%  derivatives as possible in the Taylor series beyond the deri-
%  vative of interest, in this case u1 , in order to produce a
%                                     x
%  finite difference approximation for the derivative of maximum
%  accuracy).  In this case we have only two constants, a and b,
%  to select so we can drop out only the second derivative, u1  ,
%                                                             2x
%  in the Taylor series (in addition to retaining the first deri-
%  vative).  This procedure leads to two linear algebraic equa-
%  tions in the two constants
%
%  a + 2b = 1
%
%  a + 4b = 0
%
%  Solution of these equations for a and b gives
%
%  a = 2, b = -1/2
%
%  Solution of equations (9) and (10) for u1  with these values of
%  a and b gives equation (1)               x
%
%                                      2
%  u1 = (1/2dx)(-3u1 + 4u2 - u3) + O(dx )                      (1)
%    x
%               2
%  The term O(dx ) indicates a principal error term due to trunca-
%                                                2
%  tion of the Taylor series which is of order dx .  This term in
%                    2
%  fact equals u1  dx /3f, which is easily obtained in deriving
%                3x
%  equation (1).
%
%  This same basic procedure can now be applied to the derivation
%  of equations (2) and (3).
%
%     (2)  Interior point i
%
%                                    2           3
%  a(ui-1 = ui + ui (-dx) + ui  (-dx) + ui  (-dx) + ...)
%                  x  1f      2x  2f      3x  3f
%
%                                    2           3
%  b(ui+1 = ui + ui ( dx) + ui  ( dx) + ui  ( dx) + ...)
%                  x  1f      2x  2f      3x  3f
%
%  -a + b = 1
%
%   a + b = 0
%
%  a = 1/2, b = -1/2
%                                   2
%  ui  = (1/2dx)(ui+1 - ui-1) + O(dx )                         (2)
%    x
%
%     (3)  Right end, point i = n
%
%                                      2            3
%  a(un-2 = un + un (-2dx) + un  (-2dx) + un  (-2dx) + ...)
%                  X   1f      2x   2f      3x   3f
%
%                                      2            3
%  b(un-1 = un + un ( -dx) + un  ( -dx) + un  ( -dx) + ...)
%                  x   1f      2x   2f      3x   3f
%
%  -2a - b = 1
%
%   4a + b = 0
%
%   a = -2, b = 1/2
%                                          2
%  un  = (1/2dx)(3un - 4un-1 + un-2) + O(dx )                  (3)
%    x
%
%  The weighting coefficients for equations (1), (2) and (3) can
%  be summarized as
%
%          -3   4  -1
%
%     1/2  -1   0   1
%
%           1  -4   3
%
%  Which are the coefficients reported by Bickley for n = 2, m =
%  1, p = 0, 1, 2 (Bickley, W. G., Formulae for Numerical Differ-
%  entiation, Math. Gaz., vol. 25, 1941).
%
%  Equations (1), (2) and (3) can now be programmed to generate
%  the derivative u (x) of u(x).
%                  x
%
%  Compute the spatial increment
   dx=(xu-xl)/(n-1);
   r2fdx=1./(2.*dx);
   nm1=n-1;
%
%  Equation (1) (note - the rhs of the finite difference approxi-
%  tions, equations (1), (2) and (3) have been formatted so that
%  the numerical weighting coefficients can be more easily associ-
%  ated with the Bickley matrix listed above)
   ux(1)=r2fdx*...
     (     -3.   *u(  1)     +4.   *u(  2)     -1.   *u(  3));
%
%  Equation (2)
   for i=2:nm1
     ux(i)=r2fdx*...
     (     -1.   *u(i-1)     +0.   *u(  i)     +1.   *u(i+1));
   end
%
%  Equation (3)
   ux(n)=r2fdx*...
     (      1.   *u(n-2)     -4.   *u(n-1)     +3.   *u(  n));
