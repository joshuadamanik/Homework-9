
%----------- cost/output files ---------------------------------------------

   CostDef       = 'obj_m09';                       % Cost Function Name 
 
%----------- problem definitions -------------------------------------------

   NumParamX     = 7;                             % Number of X Parameters  
   
   NumIneq       = 4;                             % Number of Inequality constraints
   NumEq     	 = 0;                             % Number of Equality constraints

   ParamMaxX = 10*[1 1 1 1 1 1 1];         % upper bounds of search space
   ParamMinX = -10*[1 1 1 1 1 1 1];		   % lower bounds of search space
 
%--------- ALM penalty & Acceptable constraint violation-----------------------------

   rho           = 10;                    % ALM penalty term  weight
   Tolerance     = 0.0001;               % constraints tolerance

%--------Annealing parameter  ----------------------------------------------------------
 
   DeciGen = 1;
   
%--------------------------------------------------------------------------

   omega = 0.995;
   phiP = 0.005;
   phiG = 0.005;