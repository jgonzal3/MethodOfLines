%  File: dss006.m
%
   function [ux]=dss006(xl,xu,n,u)
%
%  Function dss006 computes the first derivative, u , of a
%                                                  x
%  Variable u over the spatial domain xl le x le xu from classical
%  seven-point, sixth-order finite difference approximations
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
%  polynomials) are given in routines dss002 and dss004.
%
%  Seven-point formulas
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
%
%                                  2            3            4
%  f(u7 = u1 + u1 (6dx) + u1  (6dx)  + u1  (6dx)  + u1  (6dx)
%                x  1f      2x  2f       3x  3f       4x  4f
%
%                      5            6            7
%           + u1  (6dx)  + u1  (6dx)  + u1  (6dx)  + ...)
%               5x  5f       6x  6f       7x  7f
%
%  Constants a, b, c, d, e and f are selected so that the coeffi-
%  cients of the u1  terms sum to one and the coefficients of
%                  x
%  the u1  , u1  , u1  , u1   and u1   terms sum to zero.
%        2x    3x    4x    5x       6x
%
%        1      1      1      1      1
%  a +  2 b +  3 c +  4 d +  5 e +  6 f = 1
%
%        2      2      2      2      2
%  a +  2 b +  3 c +  4 d +  5 e +  6 f = 0
%
%        3      3      3      3      3
%  a +  2 b +  3 c +  4 d +  5 e +  6 f = 0
%
%        4      4      4      4      4
%  a +  2 b +  3 c +  4 d +  5 e +  6 f = 0
%
%        5      5      5      5      5
%  a +  2 b +  3 c +  4 d +  5 e +  6 f = 0
%
%        6      6      6      6      6
%  a +  2 b +  3 c +  4 d +  5 e +  6 f = 0
%
%  Simultaneous solution for a, b, c, d, e and f followed by the
%  solution of the preceding Taylor series, truncated after the
%  u   terms, for u1  gives the following seven-point approxi-
%   6x              x
%  tion
%
%  u1  = (1/6f)(-1764*u1 + 4320u2 - 5400u3 + 4800u4
%    x                                         6               (1)
%              + 2700u5 + 864u6 - 120u7) + O(dx )
%
%     (2)  Interior point, i = 2
%
%  The preceding Taylor series can be summarized for point i = 2
%  as
%
%  a(u1 = u2 + ...)
%
%  b(u3 = u2 + ...)
%
%  c(u4 = u2 + ...)
%
%  d(u5 = u2 + ...)
%
%  e(u6 = u2 + ...)
%
%  f(u7 = u2 + ...)
%
%  The corresponding algebraic equations are
%
%    1      1      1      1      1      1
%  -1 a +  1 b +  2 c +  3 d +  4 e +  5 f = 1
%
%    2      2      2      2      2      2
%  -1 a +  1 b +  2 c +  3 d +  4 e +  5 f = 0
%
%    3      3      3      3      3      3
%  -1 a +  1 b +  2 c +  3 d +  4 e +  5 f = 0
%
%    4      4      4      4      4      4
%  -1 a +  1 b +  2 c +  3 d +  4 e +  5 f = 0
%
%    5      5      5      5      5      5
%  -1 a +  1 b +  2 c +  3 d +  4 e +  5 f = 0
%
%    6      6      6      6      6      6
%  -1 a +  1 b +  2 c +  3 d +  4 e +  5 f = 0
%
%  Simultaneous solution for a, b, c, d, e and f followed by the
%  solution of the preceding Taylor series, truncated after the
%  u   terms, for u2  gives the following seven-point approxi-
%   6x              x
%  tion
%
%  u2  = (1/6f)(-120u1 - 924u2 + 1800u3 - 1200u4
%    x                                       6                 (2)
%              + 600u5 - 180u6 + 24u7) + O(dx )
%
%     (3)  Interior point, i = 3
%
%  The preceding Taylor series can be summarized for point i = 3
%  as
%
%  a(u1 = u3 + ...)
%
%  b(u2 = u3 + ...)
%
%  c(u4 = u3 + ...)
%
%  d(u5 = u3 + ...)
%
%  e(u6 = u3 + ...)
%
%  f(u7 = u3 + ...)
%
%  The corresponding algebraic equations are
%
%    1      1      1      1      1      1
%  -2 a + -1 b +  1 c +  2 d +  3 e +  4 f = 1
%
%    2      2      2      2      2      2
%  -2 a + -1 b +  1 c +  2 d +  3 e +  4 f = 0
%
%    3      3      3      3      3      3
%  -2 a + -1 b +  1 c +  2 d +  3 e +  4 f = 0
%
%    4      4      4      4      4      4
%  -2 a + -1 b +  1 c +  2 d +  3 e +  4 f = 0
%
%    5      5      5      5      5      5
%  -2 a + -1 b +  1 c +  2 d +  3 e +  4 f = 0
%
%    6      6      6      6      6      6
%  -2 a + -1 b +  1 c +  2 d +  3 e +  4 f = 0
%
%  Simultaneous solution for a, b, c, d, e and f followed by the
%  solution of the preceding Taylor series, truncated after the
%  u   terms, for u3  gives the following seven-point approxi-
%   6x              x
%  tion
%
%  u3  = (1/6f)(24u1 - 288u2 - 420u3 + 960u4
%    x                                     6                   (3)
%             - 360u5 + 96u6 - 12u7) + O(dx )
%
%     (4)  Interior point, i ne 2, 3, n-2, n-1
%
%  The preceding Taylor series can be summarized for point i = i
%  as
%
%  a(ui-3 = ui + ...)
%
%  b(ui-2 = ui + ...)
%
%  c(ui-1 = ui + ...)
%
%  d(ui+1 = ui + ...)
%
%  e(ui+2 = ui + ...)
%
%  f(ui+3 = ui + ...)
%
%  The corresponding algebraic equations are
%
%    1      1      1      1      1      1
%  -3 a + -2 b + -1 c +  1 d +  2 e +  3 f = 1
%
%    2      2      2      2      2      2
%  -3 a + -2 b + -1 c +  1 d +  2 e +  3 f = 0
%
%    3      3      3      3      3      3
%  -3 a + -2 b + -1 c +  1 d +  2 e +  3 f = 0
%
%    4      4      4      4      4      4
%  -3 a + -2 b + -1 c +  1 d +  2 e +  3 f = 0
%
%    5      5      5      5      5      5
%  -3 a + -2 b + -1 c +  1 d +  2 e +  3 f = 0
%
%    6      6      6      6      6      6
%  -3 a + -2 b + -1 c +  1 d +  2 e +  3 f = 0
%
%  Simultaneous solution for a, b, c, d, e and f followed by the
%  solution of the preceding Taylor series, truncated after the
%  u   terms, for ui  gives the following seven-point approxi-
%   6x              x
%  tion
%
%  ui  = (1/6f)(-12ui-3 + 108ui-2 - 540ui-1 + 0ui
%    x                                             6           (4)
%              + 540ui+1 - 108ui+2 + 12ui+3) + O(dx )
%
%     (5)  interior point, i = n-2
%
%  The preceding Taylor series can be summarized for point i = n-2
%  as
%
%  a(un-6 = un-2 + ...)
%
%  b(un-5 = un-2 + ...)
%
%  c(un-4 = un-2 + ...)
%
%  d(un-3 = un-2 + ...)
%
%  e(un-1 = un-2 + ...)
%
%  f(un   = un-2 + ...)
%
%  The corresponding algebraic equations are
%
%    1      1      1      1      1      1
%  -4 a + -3 b + -2 c + -1 d +  1 e +  2 f = 1
%
%    2      2      2      2      2      2
%  -4 a + -3 b + -2 c + -1 d +  1 e +  2 f = 0
%
%    3      3      3      3      3      3
%  -4 a + -3 b + -2 c + -1 d +  1 e +  2 f = 0
%
%    4      4      4      4      4      4
%  -4 a + -3 b + -2 c + -1 d +  1 e +  2 f = 0
%
%    5      5      5      5      5      5
%  -4 a + -3 b + -2 c + -1 d +  1 e +  2 f = 0
%
%    6      6      6      6      6      6
%  -4 a + -3 b + -2 c + -1 d +  1 e +  2 f = 0
%
%  Simultaneous solution for a, b, c, d, e and f followed by the
%  solution of the preceding Taylor series, truncated after the
%  u   terms, for un-2  gives the following seven-point approxi-
%   6x                x
%  tion
%
%  un-2  = (1/6f)(12un-6 - 96un-5 + 360un-4 - 960un-3
%      x                                          6            (5)
%               + 420un-2 + 288un-1 - 24un) + O(dx )
%
%     (6)  Interior point, i = n-1
%
%  The preceding Taylor series can be summarized for point i = n-1
%  as
%
%  a(un-6 = un-1 + ...)
%
%  b(un-5 = un-1 + ...)
%
%  c(un-4 = un-1 + ...)
%
%  d(un-3 = un-1 + ...)
%
%  e(un-2 = un-1 + ...)
%
%  f(un   = un-1 + ...)
%
%  The corresponding algebraic equations are
%
%    1      1      1      1      1      1
%  -5 a + -4 b + -3 c + -2 d + -1 e +  1 f = 1
%
%    2      2      2      2      2      2
%  -5 a + -4 b + -3 c + -2 d + -1 e +  1 f = 0
%
%    3      3      3      3      3      3
%  -5 a + -4 b + -3 c + -2 d + -1 e +  1 f = 0
%
%    4      4      4      4      4      4
%  -5 a + -4 b + -3 c + -2 d + -1 e +  1 f = 0
%
%    5      5      5      5      5      5
%  -5 a + -4 b + -3 c + -2 d + -1 e +  1 f = 0
%
%    6      6      6      6      6      6
%  -5 a + -4 b + -3 c + -2 d + -1 e +  1 f = 0
%
%  Simultaneous solution for a, b, c, d, e and f followed by the
%  solution of the preceding Taylor series, truncated after the
%  u   terms, for un-1  gives the following seven-point approxi-
%   6x                x
%  tion
%
%  un-1  = (1/6f)(-24un-6 + 180un-5 - 600un-4 + 1200un-3
%      x                                             6         (6)
%                - 1800un-2 + 924un-1 + 120un) + O(dx )
%
%     (7)  Right end, point i = n
%
%  The preceding Taylor series can be summarized for point i = n
%  as
%
%  a(un-6 = un   + ...)
%
%  b(un-5 = un   + ...)
%
%  c(un-4 = un   + ...)
%
%  d(un-3 = un   + ...)
%
%  e(un-2 = un   + ...)
%
%  f(un-1 = un   + ...)
%
%  The corresponding algebraic equations are
%
%    1      1      1      1      1      1
%  -6 a + -5 b + -4 c + -3 d + -2 e + -1 f = 1
%
%    2      2      2      2      2      2
%  -6 a + -5 b + -4 c + -3 d + -2 e + -1 f = 0
%
%    3      3      3      3      3      3
%  -6 a + -5 b + -4 c + -3 d + -2 e + -1 f = 0
%
%    4      4      4      4      4      4
%  -6 a + -5 b + -4 c + -3 d + -2 e + -1 f = 0
%
%    5      5      5      5      5      5
%  -6 a + -5 b + -4 c + -3 d + -2 e + -1 f = 0
%
%    6      6      6      6      6      6
%  -6 a + -5 b + -4 c + -3 d + -2 e + -1 f = 0
%
%  Simultaneous solution for a, b, c, d, e and f followed by the
%  solution of the preceding Taylor series, truncated after the
%  u   terms, for un  gives the following seven-point approxi-
%   6x              x
%  tion
%
%  un  = (1/6f)(120un-6 - 864un-5 + 2700un-4 - 4800un-3
%    x                                              6          (7)
%             + 5400un-2 - 4320un-1 + 1764un) + O(dx )
%
%  The weighting coefficients for equations (1) to (7) can be
%  summarized as
%
%            -1764   4320  -5400   4800  -2700    864   -120
%
%             -120   -924   1800  -1200    600   -180     24
%
%               24   -288   -420    960   -360     96    -12
%
%     1/6f     -12    108   -540      0    540   -108     12
%
%               12    -96    360   -960    420    288    -24
%
%              -24    180   -600   1200  -1800    924    120
%
%              120   -864   2700  -4800   5400  -4320   1764
%
%  which are the coefficients reported by Bickley for n = 6, m =
%  1, p = 0 to 6 (Bickley, W. G., Formulae for Numerical Differ-
%  entiation, Math. Gaz., vol. 25, 1941).
%
%  Equations (1) to (7) can now be programmed to generate the
%  derivative u (x) of function u(x) (arguments u and ux of
%              x
%  function dss006, respectively).
%
%  Compute the spatial increment
   dx=(xu-xl)/(n-1);
   r6fdx=1./(720.*dx);
   nm3=n-3;
%
%  Equation (1)
   ux(  1)=r6fdx*...
     ( -1764.*u(  1)  +4320.*u(  2)  -5400.*u(  3)  +4800.*u(  4)...
       -2700.*u(  5)   +864.*u(  6)   -120.*u(  7));
%
%  Equation (2)
   ux(  2)=r6fdx*...
     (  -120.*u(  1)   -924.*u(  2)  +1800.*u(  3)  -1200.*u(  4)...
        +600.*u(  5)   -180.*u(  6)    +24.*u(  7));
%
%  Equation (3)
   ux(  3)=r6fdx*...
     (   +24.*u(  1)   -288.*u(  2)   -420.*u(  3)   +960.*u(  4)...
        -360.*u(  5)    +96.*u(  6)    -12.*u(  7));
%
%  Equation (4)
   for i=4:nm3
     ux(  i)=r6fdx*...
     (   -12.*u(i-3)   +108.*u(i-2)   -540.*u(i-1)     +0.*u(  i)...
        +540.*u(i+1)   -108.*u(i+2)    +12.*u(i+3));
   end
%
%  Equation (5)
   ux(n-2)=r6fdx*...
     (   +12.*u(n-6)    -96.*u(n-5)   +360.*u(n-4)   -960.*u(n-3)...
        +420.*u(n-2)   +288.*u(n-1)    -24.*u(  n));
%
%  Equation (6)
   ux(n-1)=r6fdx*...
     (   -24.*u(n-6)   +180.*u(n-5)   -600.*u(n-4)  +1200.*u(n-3)...
       -1800.*u(n-2)   +924.*u(n-1)   +120.*u(  n));
%
%  Equation (7)
   ux(  n)=r6fdx*...
     (  +120.*u(n-6)   -864.*u(n-5)  +2700.*u(n-4)  -4800.*u(n-3)...
       +5400.*u(n-2)  -4320.*u(n-1)  +1764.*u(  n));
