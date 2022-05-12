%  File: dss008.m
%
   function [ux]=dss008(xl,xu,n,u)
%
%  Routine dss008 computes the first derivative, u , of a
%                                                 x
%  variable u over the spatial domain xl le x le xu from classical
%  nine-point, eighth-order finite difference approximations
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
%  The mathematical details of the following Taylor series (or
%  polynomials) are given in routines dss002, dss004 and
%  dss006.
%
%  Nine-point formulas
%
%     (1)  Left end, point i = 1
%
%                                  2            3            4
%  a(u2 = u1 + u1 ( dx) + u1  ( dx)  + u1  ( dx)  + u1  ( dx)
%                x  1f      2x  2f       3x  3f       4x  4f
%
%                      5            6            7
%           + u1  ( dx)  + u1  ( dx)  + u1  ( dx)  + ...)
%               5x  5f       6x  6f       7x  7f
%
%                                  2            3            4
%  b(u3 = u1 + u1 (2dx) + u1  (2dx)  + u1  (2dx)  + u1  (2dx)
%                x  1f      2x  2f       3x  3f       4x  4f
%
%                      5            6            7
%           + u1  (2dx)  + u1  (2dx)  + u1  (2dx)  + ...)
%               5x  5f       6x  6f       7x  7f
%
%                                  2            3            4
%  c(u4 = u1 + u1 (3dx) + u1  (3dx)  + u1  (3dx)  + u1  (3dx)
%                x  1f      2x  2f       3x  3f       4x  4f
%
%                      5            6            7
%           + u1  (3dx)  + u1  (3dx)  + u1  (3dx)  + ...)
%               5x  5f       6x  6f       7x  7f
%
%                                  2            3            4
%  d(u5 = u1 + u1 (4dx) + u1  (4dx)  + u1  (4dx)  + u1  (4dx)
%                x  1f      2x  2f       3x  3f       4x  4f
%
%                      5            6            7
%           + u1  (4dx)  + u1  (4dx)  + u1  (4dx)  + ...)
%               5x  5f       6x  6f       7x  7f
%
%                                  2            3            4
%  e(u6 = u1 + u1 (5dx) + u1  (5dx)  + u1  (5dx)  + u1  (5dx)
%                x  1f      2x  2f       3x  3f       4x  4f
%
%                      5            6            7
%           + u1  (5dx)  + u1  (5dx)  + u1  (5dx)  + ...)
%               5x  5f       6x  6f       7x  7f
%                                  2            3            4
%  f(u7 = u1 + u1 (6dx) + u1  (6dx)  + u1  (6dx)  + u1  (6dx)
%                x  1f      2x  2f       3x  3f       4x  4f
%
%                      5            6            7
%           + u1  (6dx)  + u1  (6dx)  + u1  (6dx)  + ...)
%               5x  5f       6x  6f       7x  7f
%
%                                  2            3            4
%  g(u8 = u1 + u1 (7dx) + u1  (7dx)  + u1  (7dx)  + u1  (7dx)
%                x  1f      2x  2f       3x  3f       4x  4f
%
%                      5            6            7
%           + u1  (7dx)  + u1  (7dx)  + u1  (7dx)  + ...)
%               5x  5f       6x  6f       7x  7f
%
%                                  2            3            4
%  h(u9 = u1 + u1 (8dx) + u1  (8dx)  + u1  (8dx)  + u1  (8dx)
%                x  1f      2x  2f       3x  3f       4x  4f
%
%                      5            6            7
%           + u1  (8dx)  + u1  (8dx)  + u1  (8dx)  + ...)
%               5x  5f       6x  6f       7x  7f
%
%  Constants a, b, c, d, e, f, g and h are selected so that the
%  coefficients of the u1  terms sum to one and the coefficients of
%                        x
%  the u1  , u1  , u1  , u1  , u1  , u1   and u1   sum to zero.
%        2x    3x    4x    5x    6x,   7x       8x
%
%        1      1      1      1      1      1      1
%  a +  2 b +  3 c +  4 d +  5 e +  6 f +  7 g +  8 h = 1
%
%        2      2      2      2      2      2      2
%  a +  2 b +  3 c +  4 d +  5 e +  6 f +  7 g +  8 h = 0
%
%        3      3      3      3      3      3      3
%  a +  2 b +  3 c +  4 d +  5 e +  6 f +  7 g +  8 h = 0
%
%        4      4      4      4      4      4      4
%  a +  2 b +  3 c +  4 d +  5 e +  6 f +  7 g +  8 h = 0
%
%        5      5      5      5      5      5      5
%  a +  2 b +  3 c +  4 d +  5 e +  6 f +  7 g +  8 h = 0
%
%        6      6      6      6      6      6      6
%  a +  2 b +  3 c +  4 d +  5 e +  6 f +  7 g +  8 h = 0
%
%        7      7      7      7      7      7      7
%  a +  2 b +  3 c +  4 d +  5 e +  6 f +  7 g +  8 h = 0
%
%        8      8      8      8      8      8      8
%  a +  2 b +  3 c +  4 d +  5 e +  6 f +  7 g +  8 h = 0
%
%  Simultaneous solution for a, b, c, d, e, f, g and h followed by
%  the solution of the preceding Taylor series, truncated after the
%  u   terms, for u1  gives the following nine-point approximation
%   8x              x
%
%  u1  = (1/8f)(-109584u1 + 322560u2 - 564480u3 + 752640u4
%    x
%               -705600u5 + 451584u6 - 188160u7 + 46080u8
%                               8
%               - 5040u9) + O(dx )                             (1)
%
%  The preceding analysis can be repeated to produce nine-point
%  approximations for the first derivatives u2 , u3 , u4 , ui ,
%                                             x    x    x    x
%  un-3 , un-2 , un-1  and un .  The results can be summarized by
%      x      x      x       x
%  the following Bickley matrix for n = 8, m = 1, p = 0 to 8,
%  (Bickley, W. G., Formulae for Numerical Differentiation, Math.
%  Gaz., vol. 25, 1941)
%
%            -109584   322560  -564480   752640  -705600
%
%              -5040   -64224   141120  -141120   117600
%
%                720   -11520   -38304    80640   -50400
%
%               -240     2880   -20160   -18144    50400
%
%     1/8f       144    -1536     8064   -32256        0
%
%               -144     1440    -6720    20160   -50400
%
%                240    -2304    10080   -26880    50400
%
%               -720     6720   -28224    70560  -117600
%
%               5040   -46080   188160  -451584   705600
%
%  From this Bickley matrix, the finite difference approximation of
%  the first derivative can be programmed for each of the grid points
%  1, 2, 3, 4,..., i,..., n-3, n-2, n-1, n (taking into account
%  the symmetry properties of the Bickley matrix).
%
%  Compute the spatial increment
   dx=(xu-xl)/(n-1);
   r8fdx=1./(40320.*dx);
   nm4=n-4;
%
%  Grid point 1
   ux(  1)=r8fdx*...
     (-109584.     *u(  1)...
      +322560.     *u(  2)...
      -564480.     *u(  3)...
      +752640.     *u(  4)...
      -705600.     *u(  5)...
      +451584.     *u(  6)...
      -188160.     *u(  7)...
      +46080.      *u(  8)...
      -5040.       *u(  9));
%
%  Grid point 2
   ux(  2)=r8fdx*...
     (-5040.       *u(  1)...
      -64224.      *u(  2)...
      +141120.     *u(  3)...
      -141120.     *u(  4)...
      +117600.     *u(  5)...
      -70560.      *u(  6)...
      +28224.      *u(  7)...
      -6720.       *u(  8)...
      +720.        *u(  9));
%
%  Grid point 3
   ux(  3)=r8fdx*...
     (+720.        *u(  1)...
      -11520.      *u(  2)...
      -38304.      *u(  3)...
      +80640.      *u(  4)...
      -50400.      *u(  5)...
      +26880.      *u(  6)...
      -10080.      *u(  7)...
      +2304.       *u(  8)...
      -240.        *u(  9));
%
%  Grid point 4
   ux(  4)=r8fdx*...
     (-240.        *u(  1)...
      +2880.       *u(  2)...
      -20160.      *u(  3)...
      -18144.      *u(  4)...
      +50400.      *u(  5)...
      -20160.      *u(  6)...
      +6720.       *u(  7)...
      -1440.       *u(  8)...
      +144.        *u(  9));
%
%  Grid point i
   for i=5:nm4
     ux(  i)=r8fdx*...
     (+144.        *u(i-4)...
      -1536.       *u(i-3)...
      +8064.       *u(i-2)...
      -32256.      *u(i-1)...
      +0.          *u(i  )...
      +32256.      *u(i+1)...
      -8064.       *u(i+2)...
      +1536.       *u(i+3)...
      -144.        *u(i+4));
   end
%
%  Grid point n-3
   ux(n-3)=r8fdx*...
     (-144.        *u(n-8)...
      +1440.       *u(n-7)...
      -6720.       *u(n-6)...
      +20160.      *u(n-5)...
      -50400.      *u(n-4)...
      +18144.      *u(n-3)...
      +20160.      *u(n-2)...
      -2880.       *u(n-1)...
      +240.        *u(n  ));
%
%  Grid point n-2
   ux(n-2)=r8fdx*...
     (+240.        *u(n-8)...
      -2304.       *u(n-7)...
      +10080.      *u(n-6)...
      -26880.      *u(n-5)...
      +50400.      *u(n-4)...
      -80640.      *u(n-3)...
      +38304.      *u(n-2)...
      +11520.      *u(n-1)...
      -720.        *u(n  ));
%
%  Grid point n-1
   ux(n-1)=r8fdx*...
     (-720.        *u(n-8)...
      +6720.       *u(n-7)...
      -28224.      *u(n-6)...
      +70560.      *u(n-5)...
      -117600.     *u(n-4)...
      +141120.     *u(n-3)...
      -141120.     *u(n-2)...
      +64224.      *u(n-1)...
      +5040.       *u(n  ));
%
%  Grid point n
   ux(n  )=r8fdx*...
     (+5040.       *u(n-8)...
      -46080.      *u(n-7)...
      +188160.     *u(n-6)...
      -451584.     *u(n-5)...
      +705600.     *u(n-4)...
      -752640.     *u(n-3)...
      +564480.     *u(n-2)...
      -322560.     *u(n-1)...
      +109584.     *u(n  ));
