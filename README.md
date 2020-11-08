## File Descriptions

#### dogleg.m
To run this function call dogleg(xn, fx, gx, H, H_, delta)\
[although is called as part of tr_dogleg]\\
where
- xn : initial point
- fx : actual function
- gx : gradient function
- H  : Hessian matrix
- H_ : inverse Hessian matrix
- delta : radius of the trust region

#### sr1.m
To run this function call sr1(H, H_, d, y, eta)\
[although is called as part of tr_dogleg]\\
where
- H  : Hessian matrix
- H_ : inverse Hessian matrix
- d  : step between two points
- y  : difference between gradient values
- eta : small value greater than zero

#### tr_dogleg.m
To run this function call tr_dogleg(f,df,H,H_,xn,delta,delmax,rho_ac,tol)\
where
- f  : actual function
- df : gradient function
- H  : Hessian matrix
- H_ : inverse Hessian matrix
- xn : initial point
- delta : radius of the trust region
- delmax : maximum radius of the trust region
- rho_ac : accuracy
- tol : tolerance

Note that eta, is defined as 10^(-6) in the file tr_dogleg.m
