
%----------- cost/output files ---------------------------------------------

   CostDef       = 'obj_m10';                       % Cost Function Name 

%----------- problem definitions -------------------------------------------

   NumParamX     = 8;                            % Number of X Parameters  
 
   NumIneq       = 6;                             % Number of Inequality constraints
   NumEq     	 = 0;                             % Number of Equality constraints

   ParamMaxX = 10000*[1 1 1 1 1 1 1 1];            % upper bounds of search space
   ParamMinX = 10*[10 100 100 1 1 1 1 1];		   % lower bounds of search space

%--------- ALM penalty & Acceptable constraint violation-----------------------------

   rho           = 10;                    % ALM penalty term  weight
   Tolerance     = 0.0001;               % constraints tolerance

%--------Annealing parameter  ----------------------------------------------------------
 
   DeciGen = 3000;

% recommended:  DeciGen > 3000, MaxGen > 15000   
%--------------------------------------------------------------------------

   omega = 0.995;
   phiP = 0.005;
   phiG = 0.005;