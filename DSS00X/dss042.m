%  File: dss042.m
%
   function [uxx]=dss042(xl,xu,n,u,ux,nl,nu)
%
%  Function dss042 computes a second-order approximation of a
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
%             x = xu (input).  the allowable values are
%
%                1 - Dirichlet boundary condition at x = xu
%                    (ux(n) is not used)
%
%                2 - Neumann boundary condition at x = xu
%                    (ux(n) is used)
%
%  The following derivation is for a set of second-order, four-point
%  approximations for a second derivative that can be used at the
%  boundaries of a spatial domain.  These approximations have the
%  features
%
%     (1)  Only interior and boundary points are used (i.e., no
%          fictitious points are used)
%
%     (2)  The normal derivative at the boundary is included as part
%          of the approximation for the second derivative
%
%     (3)  Approximations for the boundary conditions are not used.
%
%  The derivation is by Professor Gilbert A. Stengle, Department of
%  Mathematics, Lehigh University, Bethlehem, pa 18015, and was done
%  on December 7, 1985.
%
%  For an approximation at the left boundary, involving the points
%  i, i+1, i+2 and i+3, consider the following Taylor series expan-
%  sions
%
%                  ux(i)( dx)   uxx(i)( dx)**2   uxxx(i)( dx)**3
%  u(i+1) = u(i) + ---------- + -------------- + --------------- +...
%                       1             2                 6
%
%
%                  ux(i)(2dx)   uxx(i)(2dx)**2   uxxx(i)(2dx)**3
%  u(i+2) = u(i) + ---------- + -------------- + --------------- +...
%                       1             2                 6
%
%  If we now form the following linear combination, involving con-
%  stants a, b, c and d to be determined, and use the preceding two
%  Taylor series,
%
%     a*u(i) + b*ux(i) + c*u(i+1) + d*u(i+2)
%
%  we have
%
%     a*u(i) + b*ux(i) + c*u(i+1) + d*u(i+2) =
%
%     (a + b + c + d)*u(i) +
%
%     (b + dx*c + 2*dx*d)*ux(i) +
%
%     (c*(dx**2)/2 + d*((2*dx)**2)/2)*uxx(i) +
%
%     (c*(dx**3)/6 + d*((2*dx)**3)/6)*uxxx(i) + O(dx**4)
%
%  The third derivative, uxxx(i), can be dropped by taking
%
%     c = -8*d
%
%  The second derivative, uxx(i), can be retained by taking
%
%     (dx**2)(c/2 + 2*d) = 1
%
%  which, when combined with the preceding result gives
%
%     d = -1/(2*(dx**2))
%
%     c = 4/(dx**2)
%
%  The first derivative, ux(i), can be dropped by taking
%
%     b + dx*c + 2*dx*d = 0
%
%  or
%
%     b = -dx*c - 2*dx*d = -4/dx - 2*dx*(-1/(2*(dx**2))) = -3/dx
%
%  Finally, u(i), can be dropped by taking
%
%     a = - c - d = 8*d - d = -7*d = -7/(2*(dx**2))
%
%  If we now solve for the derivative of interest, uxx(i),
%
%     uxx(i) = -7/(2(dx**2))*u(i) - 3/dx*ux(i)
%
%              + 8/(dx**2)*u(i+1) - 1/(2*(dx**2))u(i+2) + O(dx**2)
%
%       = (1/(2*(dx**2)))*(-u(i+2) + 8*u(i+1) - 7*u(i) - 6*dx*ux(i))
%
%         + O(dx**2)
%
%  which is the four-point, second-order approximation for the second
%  derivative, uxx(i), including the first derivative, ux(i).
%
%  Four checks of this approximation can easily be made for u(i) =
%  1, u(i) = x, u(i) = x**2 and u(i) = x**3
%
%     uxx(i) = (1/(2*(dx**2)))*(-1 + 8*1 - 7*1 - 6*dx*0) = 0
%
%     uxx(i) = (1/(2*(dx**2)))*(-(x + 2*dx) + 8*(x + dx)
%
%              -7*x - 6*dx*1) = 0
%
%     uxx(i) = (1/(2*(dx**2)))*(-(x + 2*dx)**2 + 8*(x + dx)**2
%
%            - 7*(x**2) - 6*dx*(2*x))
%
%             = (-  x**2 -  4*x*dx - 4*dx**2
%
%               + 8*x**2 + 16*x*dx + 8*dx**2
%
%               - 7*x**2 - 12*x*dx)/(2*(dx**2)) = 2
%
%     uxx(i) = (1/(2*(dx**2)))*(-(x + 2*dx)**3 + 8*(x + dx)**3
%
%            - 7*(x**3) - 6*dx*(3*x**2))
%
%            = (1/(2*(dx**2)))*(- x**3 - 6*dx*x**2 - 12*x*dx**2
%
%            - 8*dx**3 + 8*x**3 + 24*dx*x**2 + 24*x*dx**2 + 8*dx**3
%
%            - 7*x**3 - 18*dx*x**2)
%
%            = (1/(2*(dx**2)))*(12*x*dx**2) = 6*x
%
%  The preceding approximation for uxx(i) can be applied at the
%  left boundary value of x by taking i = 1.  An approximation at
%  the right boundary is obtained by taking dx = -dx and reversing
%  the subscripts in the preceding approximation, with i = n
%
%     uxx(i)
%
%       = (1/(2*(dx**2)))*(-u(i-2) + 8*u(i-1) - 7*u(i) + 6*dx*ux(i))
%
%         + O(dx**2)
%
%  To obtain approximations of the second derivative which do not
%  involve the first derivative, we take as the linear combination
%
%     a*u(i) + b*u(i+1) + c*u(i+2) + d*u(i+3) 
%
%  we have
%
%     a*u(i) + b*u(i+1) + c*u(i+2) + d*u(i+3) =
%
%     (a + b + c + d)*u(i)+
%
%     (dx*b + 2*dx*c + 4*dx*d)*ux(i)+
%
%     (b*(dx**2)/2 + c*((2*dx)**2)/2 + d*((3*dx)**2)/2)*uxx(i) +
%
%     (b*(dx**3)/6 + c*((2*dx)**3)/6 + d*((3*dx)**3)/6)*uxx(i) +
%
%     O(dx**4)
%
%  The third derivative, uxxx(i), can be dropped by taking
%
%     b + 8*c + 27*d = 0
%
%  The second derivative, uxx(i), can be retained by taking
%
%     (dx**2)*(b/2 + 2*c + (9/2)*d) = 1
%
%  The first derivative can be dropped by taking
%
%     b + 2*c + 3*d = 0
%
%  Solution of the preceding equations for c and d by elimination of
%  b gives
%
%     6*c + 24*d = 0
%
%     4*c + 18*d = -2/(dx**2)
%
%  Then, eliminating c, gives
%
%     (18 - 16)*d = -2/(dx**2)
%
%  or
%
%     d = -1/(dx**2)
%
%     c = (24/6)/(dx**2) = 4/(dx**2)
%
%     b = -8/(dx**2) + 3/(dx**2) = -5/(dx**2)
%
%  u(i) can be dropped by taking
%
%     a + b + c + d = 0
%
%  or
%
%     a = (5 - 4 + 1)/(dx**2) = 2/(dx**2)
%
%  If we now solve for the derivative of interest, uxx(i),
%
%     uxx(i) = (1/dx**2)*(2*u(i) - 5*u(i+1) + 4*u(i+2) - 1*u(i+3))
%
%            + O(dx**2)
%
%  Which is the four-point, second-order approximation for the second
%  derivative, uxx(i), without the first derivative, ux(i).
%
%  Four checks of this approximation can easily be made for u(i) =
%  1, u(i) = x, u(i) = x**2 and u(i) = x**3
%
%     uxx(i) = (1/dx**2)*(2 - 5 + 4 - 1) = 0
%
%     uxx(i) = (1/dx**2)*(2*x - 5*(x + dx) + 4*(x + 2*dx)
%
%              - 1*(x + 3*dx)) = 0
%
%     uxx(i) = (1/dx**2)*(2*x**2 - 5*(x + dx)**2 + 4*(x + 2*dx)**2
%
%              - 1*(x + 3*dx)**2) = 2
%
%     uxx(i) = (1/dx**2)*(2*x**3 - 5*(x + dx)**3 + 4*(x + 2*dx)**3
%
%            - 1*(x + 3*dx)**3)
%
%             = (1/dx**2)*(2*x**3 - 5*x**3 - 15*x*dx**2
%
%             - 15*dx*x**2 - 5*dx**3 + 4*x**3 + 24*dx*x**2
%
%             + 48*x*dx**2 + 32*dx**3 - x**3 - 9*dx*x**2
%
%             - 27*x*dx**2 - 27dx**3)
%
%             = (1/dx**2)*(6*x*dx**2) = 6*x
%
%  The preceding approximation for uxx(i) can be applied at the
%  left boundary value of x by taking i = 1.  An approximation at
%  the right boundary is obtained by taking dx = -dx and reversing
%  the subscripts in the preceding approximation, with i = n
%
%     uxx(i) = (1/dx**2)*(2*u(i) - 5*u(i-1) + 4*u(i-2) - 1*u(i-3))
%
%            + O(dx**2)
%
%  Grid spacing
   dx=(xu-xl)/(n-1);
%
%  uxx at the left boundary, without ux
   if nl==1
     uxx(1)=(( 2.)*u(  1)...
            +(-5.)*u(  2)...
            +( 4.)*u(  3)...
            +(-1.)*u(  4))/(dx^2);
%
%  uxx at the left boundary, including ux
   elseif nl==2
     uxx(1)=((-7.)*u(  1)...
            +( 8.)*u(  2)...
            +(-1.)*u(  3))/(2.*dx^2)...
            +(-6.)*ux( 1) /(2.*dx);
      end
%
%  uxx at the right boundary, without ux
   if nu==1
     uxx(n)=(( 2.)*u(n  )...
            +(-5.)*u(n-1)...
            +( 4.)*u(n-2)...
            +(-1.)*u(n-3))/(dx^2);
%
%  uxx at the right boundary, including ux
   elseif nu==2
     uxx(n)=((-7.)*u(n  )...
            +( 8.)*u(n-1)...
            +(-1.)*u(n-2))/(2.*dx^2)...
            +( 6.)*ux(n ) /(2.*dx);
   end
%
%  uxx at the interior grid points
   for i=2:n-1
     uxx(i)=(u(i+1)-2.*u(i)+u(i-1))/dx^2;
   end
