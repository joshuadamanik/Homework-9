
%----------- cost/output files ---------------------------------------------

   CostDef       = 'obj_m01';                       % Cost Function Name 

%----------- problem definitions -------------------------------------------

   NumParamX     = 13;                            % Number of X Parameters  

   NumIneq       = 9;                             % Number of Inequality constraints
   NumEq     	 = 0;                             % Number of Equality constraints

   ParamMaxX = [1 1 1 1 1 1 1 1 1 100 100 100 1];       % upper bounds of search space
   ParamMinX = [0 0 0 0 0 0 0 0 0 0   0    0  0];	    % lower bounds of search space
 
%--------- ALM penalty & Acceptable constraint violation-----------------------------

   rho           = 10;                    % ALM penalty term  weight
   Tolerance     = 0.0001;               % constraints tolerance

%--------Annealing parameter  ----------------------------------------------------------
 
   DeciGen = 1;

%--------------------------------------------------------------------------

   omega = 0.995;
   phiP = 0.005;
   phiG = 0.005;