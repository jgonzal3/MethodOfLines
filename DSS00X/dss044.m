%  File: ds044.m
%
   function [uxx]=dss044(xl,xu,n,u,ux,nl,nu)
%
%  Function dss044 computes a fourth-order approximation of a
%  second-order derivative, with or without the normal derivative
%  at the boundary.
%
%  Argument list
%
%     xl      Left value of the spatial independent variable (input)
%
%     xu      Right value of the spatial independent variable (input)
%
%     n       Number of spatial grid points, including the end
%             points (input)
%
%     u       One-dimensional array of the dependent variable to be
%             differentiated (input)
%
%     ux      One-dimensional array of the first derivative of u.
%             The end values of ux, ux(1) and ux(n), are used in
%             Neumann boundary conditions at x = xl and x = xu,
%             depending on the arguments nl and nu (see the de-
%             scription of nl and nu below)
%
%     uxx     One-dimensional array of the second derivative of u
%             (output)
%
%     nl      Integer index for the type of boundary condition at
%             x = xl (input).  The allowable values are
%
%                1 - Dirichlet boundary condition at x = xl
%                    (ux(1) is not used)
%
%                2 - Neumann boundary condition at x = xl
%                    (ux(1) is used)
%
%     nu      Integer index for the type of boundary condition at
%             x = xu (input).  The allowable values are
%
%                1 - Dirichlet boundary condition at x = xu
%                    (ux(n) is not used)
%
%                2 - Neumann boundary condition at x = xu
%                    (ux(n) is used)
%
%  The following derivation was completed by W. E. Schiesser, Depts
%  of CHE and Math, Lehigh University, Bethlehem, PA 18015, USA, on
%  December 15, 1986.  Additional details are given in function
%  dss042.
%
%  ******************************************************************
%
%  (1)  uxx at the interior points 3, 4,..., n-2
%
%  To develop a set of fourth-order correct differentiation formulas
%  for the second derivative uxx, we consider first the interior
%  grid points at which a symmetric formula can be used.
%
%  If we consider a formula of the form
%
%     a*u(i-2) + b*u(i-1) + e*u(i) + c*u(i+1) + d*u(i+2)
%
%  Taylor series expansions of u(i-2), u(i-1), u(i+1) and u(i+2)
%  can be substituted into this formula.  We then consider the
%  linear albegraic equations relating a, b, c and d which will
%  retain certain terms, i.e., uxx, and drop others, e.g., uxxx,
%  uxxxx and uxxxxx.
%
%  Thus, for grid points 3, 4,..., n-2
%
%     To retain uxx
%
%        4*a + b + c +  4*d = 2                              (1)
%
%     To drop uxxx
%
%       -8*a - b + c +  8*d = 0                              (2)
%
%     To drop uxxxx
%
%       16*a + b + c + 16*d = 0                              (3)
%
%     To drop uxxxxx
%
%      -32*a - b + c + 32*d = 0                              (4)
%
%  Equations (1) to (4) can be solved for a, b, c and d.  If equa-
%  tion (1) is added to equation (2)
%
%        -4*a + 2*c + 12*d = 2                               (5)
%
%  If equation (1) is subtracted from equation (3)
%
%        12*a + 12*d = -2                                    (6)
%
%  If equation (1) is added to equation (4)
%
%        -28*a + 2*c + 36*d = 2                              (7)
%
%  Equations (5) to (7) can be solved for a, c and d.  If equation
%  (5) is subtracted from equation (7), and the result combined
%  with equation (6)
%
%         12*a + 12*d = -2                                   (6)
%
%        -24*a + 24*d = 0                                    (8)
%
%  Equations (6) and (8) can be solved for a and d.  From (8), a
%  = d.  From equation (6), a = -1/12 and d = -1/12.  Then, from
%  equation (5), c = 4/3, and from equation (1), b = 4/3.
%
%  The final differentiation formula is then obtained as
%
%     (-1/12)*u(i-2) +   (4/3)*u(i-1) +
%
%       (4/3)*u(i+1) + (-1/12)*u(i+2)
%
%     (-1/12 + 4/3 - 1/12 + 4/3)*u(i) + uxx(i)*(dx**2) + O(dx**6)
%
%  or
%
%     uxx(i) = (1/(12*dx**2))*(-1*u(i-2) + 16*u(i-1)
%
%                  - 30*u(i) + 16*u(i+1) -  1*u(i+2)         (9)
%
%                  + O(dx**4)
%
%  Note that the ux term drops out, i.e., the basic equation is
%
%        -2*a - b + c + 2*d =
%
%        -2*(-1/12) - (4/3) + (4/3) + 2*(-1/12) = 0
%
%  Equation (9) was obtained by dropping all terms in the underlying
%  Taylor series up to and including the fifth derivative, uxxxxx.
%  Thus, equation (9) is exact for polynomials up to and including
%  fifth order.  This can be checked by substituting the functions
%  1, x, x**2, x**3, x**4 and x**5 in equation (9) and computing the
%  corresponding derivatives for comparison with the known second
%  derivatives.  This is done for 1 merely by summing the weighting
%  coefficients in equation (9), which should sum to zero, i.e.,
%  -1 + 16 - 30 + 16 -1 = 0.
%
%  For the remaining functions, the algebra is rather involved, but
%  these functions can be checked numerically, i.e., numerical values
%  of x**2, x**3, x**4 and x**5 can be substituted in equation (9)
%  and the computed derivatives can be compared with the know numeri-
%  cal second derivatives.  This is not a proof of correctness of
%  equation (9), but would likely detect any errors in equation (9). 
%
%  ******************************************************************
%
%  (2)  uxx at the interior points i = 2 and n-1
%
%  For grid point 2, we consider a formula of the form
%
%     a*u(i-1) + f*u(i) + b*u(i+1) + c*u(i+2) + d*u(i+3) + e*u(i+4)
%
%  Taylor series expansions of u(i-1), u(i+1), u(i+2), u(i+3) and
%  u(i+4) when substituted into this formula give linear algebraic
%  equations relating a, b, c, d and e.
%
%     To drop ux
%
%        -a + b + 2*c + 3*d + 4*e = 0                       (10)
%
%     To retain uxx
%
%        a + b + 4*c + 9*d + 16*e = 2                       (11)
%
%     To drop uxxx
%
%       -a + b + 8*c + 27*d + 64*e = 0                      (12)
%
%     To drop uxxxx
%
%        a + b + 16*c + 81*d + 256*e = 0                    (13)
%
%     To drop uxxxxx
%
%       -a + b + 32*c + 243*d + 1024*e = 0                  (14)
%
%  Equations (11), (12), (13) and (14) can be solved for a, b, c,
%  d and e.  If equation (10) is added to equation (11)
%
%     2*b + 6*c + 12*d +20*e = 2                            (15)
%
%  If equation (10) is subtracted from equation (12)
%
%     6*c + 24*d + 60*e = 0                                 (16)
%
%  If equation (10) is added to equation (13)
%
%     2*b + 18*c + 84*d + 260*e = 0                         (17)
%
%  If equation (10) is subtracted from equation (14)
%
%     30*c + 240*d + 1020*e = 0                             (18)
%
%  Equations (15), (16), (17) and (18) can be solved for b, c, d
%  and e.
%
%     6*c + 24*d + 60*e = 0                                 (16)
%
%  If equation (15) is subtracted from equation (17)
%
%     12*c + 72*d + 240*e = -2                              (19)
%
%     30*c + 240*d + 1020*e = 0                             (18)
%
%  Equations (16), (18) and (19) can be solved for c, d and e.  If
%  two times equation (16) is subtracted from equation (19),
%
%     24*d + 120*e = -2                                     (20)
%
%  If five times equation (16) is subtracted from equation (18),
%
%     120*d + 720*e = 0                                     (21)
%
%  Equations (20) and (21) can be solved for d and e.  From (21),
%  e = (-1/6)*d.  Substitution in equation (20) gives d = -1/2.
%  thus, e = 1/12.  From equation (16), c = 7/6.  From equation
%  (15), b = -1/3.  From equation (10), a = 5/6.
%
%  The final differentiation formula is then obtained as
%
%  (5/6)*u(i-1) + (-1/3)*u(i+1) + (7/6)*u(i+2) + (-1/2)*u(i+3)
%
%  + (1/12)*u(i+4) = (5/6 - 1/3 + 7/6 - 1/2 + 1/12)*u(i)
%
%  + uxx*(dx**2) + O(dx**6)
%
%  or
%
%  uxx(i) = (1/12*dx**2)*(10*u(i-1) - 15*u(i) - 4*u(i+1)
%                                                           (22)
%         + 14*u(i+2) - 6*u(i+3) + 1*u(i+4)) + O(dx**4)
%
%  Equation (22) will be applied at i = 2 and n-1.  thus
%
%  uxx(2) = (1/12*dx**2)*(10*u(1) - 15*u(2) - 4*u(3)
%                                                           (23)
%         + 14*u(4) - 6*u(5) + 1*u(6)) + O(dx**4)
%
%  uxx(n-1) = (1/12*dx**2)*(10*u(n) - 15*u(n-1) - 4*u(n-2)
%                                                           (24)
%           + 14*u(n-3) - 6*u(n-4) + 1*u(n-5)) + O(dx**4) 
%
%  ******************************************************************
%
%  (3)  uxx at the boundary points 1 and n
%
%  Finally, for grid point 1, an approximation with a Neumann bound-
%  ary condition of the form
%
%     a*u(i+1) + b*u(i+2) + c*u(i+3) + d*u(i+4) + e*ux(i) + f*u(i)
%
%  Will be used.  the corresponding algebraic equations are
%
%     To drop ux
%
%        a + 2*b + 3*c + 4*d + e = 0                        (25)
%
%     To retain uxx
%
%        a + 4*b + 9*c + 16*d = 2                           (26)
%
%     To drop uxxx
%
%        a + 8*b + 27*c + 64*d = 0                          (27)
%
%     To drop uxxxx
%
%        a + 16*b + 81*c + 256*d = 0                        (28)
%
%     To drop uxxxxx
%
%        a + 32*b + 243*c + 1024*d = 0                      (29)
%
%  Equations (25) to (29) can be solved for a, b, c, d and e.  If
%
%  Equation (26) is subtracted from equations (27), (28) and (29),
%
%     4*b + 18*c + 48*d = -2                                (30)
%
%     12*b + 72*c + 240*d = -2                              (31)
%
%     28*b + 234*c + 1008*d = -2                            (32)
%
%  Equations (30), (31) and (32) can be solved for b, c and d
%
%     18*c + 96*d = 4                                       (33)
%
%     108*c + 672*d = 12                                    (34)
%
%  Equations (3) and (34) can be solved for c and d, c = 8/9, d =
%  -1/8.
%
%  From equation (30), b = -3.  From equation (26), a = 8.  From 
%  equation (25), e = -25/6.
%
%  The final differentiation formula is then obtained as
%
%  8*u(i+1) - 3*u(i+2) + (8/9)*u(i+3) - (1/8)*u(i+4)
%
%  - (25/6)*ux(i)*dx
%
%  = (8 - 3 + (8/9) - (1/8))*u(i) + uxx*(dx**2) + O(dx**6)
%
%  or
%
%  uxx(i) = (1/12*dx**2)*((-415/6)*u(i) + 96*u(i+1) - 36*u(i+2)
%                                                                (35)
%  + (32/3)*u(i+3) - (3/2)*u(i+4) - 50*ux(i)*dx) + O(dx**4)
%
%  Equation (35) will be applied at i = 1 and i = n
%
%  uxx(1) = (1/12*dx**2)*((-415/6)*u(1) + 96*u(2) - 36*u(3)
%                                                                (36)
%  + (32/3)*u(4) - (3/2)*u(5) - 50*ux(1)*dx) + O(dx**4)
%
%  uxx(n) = (1/12*dx**2)*((-415/6)*u(n) + 96*u(n-1) - 36*u(n-2)
%                                                                (37)
%  + (32/3)*u(n-3) - (3/2)*u(n-4) + 50*ux(n)*dx) + O(dx**4) 
%
%  Alternatively, for grid point 1, an approximation with a Dirichlet
%  boundary condition of the form
%
%  a*u(i+1) + b*u(i+2) + c*u(i+3) + d*u(i+4) + e*u(i+5) + f*u(i)
%
%  can be used.  The corresponding algebraic equations are
%
%     To drop ux
%
%        a + 2*b + 3*c + 4*d + 5*e = 0                      (38)
%
%     To retain uxx
%
%        a + 4*b + 9*c + 16*d + 25*e = 2                    (39)
%
%     To drop uxxx
%
%        a + 8*b + 27*c + 64*d + 125*e = 0                  (40)
%
%     To drop uxxxx
%
%        a + 16*b + 81*c + 256*d + 625*e = 0                (41)
%
%     To drop uxxxxx
%
%        a + 32*b + 243*c + 1024*d + 3125*e = 0             (42)
%
%  Equations (38), (39), (40), (41) amd (42) can be solved for a,
%  b, c, d and e.
%
%        2*b + 6*c + 12*d + 20*e = 2                        (43)
%
%        6*b + 24*c + 60*d + 120*e = 0                      (44)
%
%        14*b + 78*c + 252*d + 620*e = 0                    (45)
%
%        30*b + 240*c + 1020*d + 3120*e = 0                 (46)
%
%  Equations (43), (44), (45) and (46) can be solved for b, c, d
%  and e
%
%        6*c + 24*d + 60*e = -6                             (47)
%
%        36*c + 168*d + 480*e = -14                         (48)
%
%        150*c + 840*d + 2820*e = -30                       (49)
%
%  Equations (47), (48) and (49) can be solved for c, d and e
%
%        24*d + 120*e = 22                                  (50)
%
%        240*d + 1320*e = 120                               (51)
%
%  From equations (50) and (51), d = 61/12, e = -5/6.  From equation 
%  (47), c = -13.  From equation (43), b = 107/6.  From equation (38), 
%  a = -77/6.
%
%  The final differentiation formula is then obtained as
%
%  (-77/6)*u(i+1) + (107/6)*u(i+2) - 13*u(i+3) + (61/12)*u(i+4)
%
%  - (5/6)*u(i+5) = (-77/6 + 107/6 - 13 + 61/12 - 5/6)*u(i) +
%
%  uxx(i)*(dx**2) + O(dx**6)
%
%  or
%
%  uxx(i) = (1/12*dx**2)*(45*u(i) - 154*u(i+1) + 214*u(i+2)
%                                                                (52)
%         - 156*u(i+3) + 61*u(i+4) - 10*u(i+5)) + O(dx**4)
%
%  Equation (52) will be applied at i = 1 and i = n
%
%  uxx(1) = (1/12*dx**2)*(45*u(1) - 154*u(2) + 214*u(3)
%                                                                (53)
%         - 156*u(4) + 61*u(5) - 10*u(6)) + O(dx**4)
%
%  uxx(n) = (1/12*dx**2)*(45*u(n) - 154*u(n-1) + 214*u(n-2)
%                                                                (54)
%         -156*u(n-3) + 61*u(n-4) - 10*u(n-5)) + O(dx**4) 
%
%  ******************************************************************
%
%  Grid spacing
   dx=(xu-xl)/(n-1);
%
%  1/(12*dx**2) for subsequent use
   r12dxs=1./(12.0*dx^2);
%
%  uxx at the left boundary
%
%     Without ux (equation (53))
      if nl==1
        uxx(1)=r12dxs*...
                      (      45.0*u(1)...
                           -154.0*u(2)...
                           +214.0*u(3)...
                           -156.0*u(4)...
                            +61.0*u(5)...
                            -10.0*u(6));
%
%     With ux (equation (36))
      elseif nl==2
        uxx(1)=r12dxs*...
                      (-415.0/6.0*u(1)...
                            +96.0*u(2)...
                            -36.0*u(3)...
                        +32.0/3.0*u(4)...
                         -3.0/2.0*u(5)...
                              -50.0*ux(1)*dx);
         end
%
%  uxx at the right boundary
%
%     Without ux (equation (54))
      if nu==1
        uxx(n)=r12dxs*...
                      (      45.0*u(n  )...
                           -154.0*u(n-1)...
                           +214.0*u(n-2)...
                           -156.0*u(n-3)...
                            +61.0*u(n-4)...
                            -10.0*u(n-5));
%
%     With ux (equation (37))
      elseif nu==2
        uxx(n)=r12dxs*...
                      (-415.0/6.0*u(n  )...
                            +96.0*u(n-1)...
                            -36.0*u(n-2)...
                        +32.0/3.0*u(n-3)...
                         -3.0/2.0*u(n-4)...
                            +50.0*ux(n  )*dx);
      end
%
%  uxx at the interior grid points
%
%     i = 2 (equation (23))
         uxx(2)=r12dxs*...
                       (     10.0*u(1)...
                            -15.0*u(2)...
                             -4.0*u(3)...
                            +14.0*u(4)...
                             -6.0*u(5)...
                             +1.0*u(6));
%
%     i = n-1 (equation (24))
         uxx(n-1)=r12dxs*...
                         (   10.0*u(n  )...
                            -15.0*u(n-1)...
                             -4.0*u(n-2)...
                            +14.0*u(n-3)...
                             -6.0*u(n-4)...
                             +1.0*u(n-5));
%
%     i = 3, 4,..., n-2 (equation (9))
         for i=3:n-2
         uxx(i)=r12dxs*...
                       (     -1.0*u(i-2)...
                            +16.0*u(i-1)...
                            -30.0*u(i  )...
                            +16.0*u(i+1)...
                             -1.0*u(i+2));
         end
