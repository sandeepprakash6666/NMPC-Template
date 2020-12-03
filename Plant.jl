
using DifferentialEquations
using Sundials
using Plots

include("parameters.jl")

## Model Equations for Plant (Change here)
function Plant_Eqns(x,z,u)

		#region-> Debugging (test Arguments)
			# x = [0.0; 0.0]
			# z = [0.53]
			# u = [0.0]
		#endregion

		#region-> Loading Parameters
		model_par = Model_par()
			x2_f  = model_par[:x2_f]
			km    = model_par[:km]
			k1    = model_par[:k1]
			Y     = model_par[:Y]
			μ_max = model_par[:μ_max]

			dx  = 0*copy(x)    #this creates a copy, not a pointer
			alg = 0*copy(z)
		#endregion

		#region->renaming variables to write equations more naturally
			x1 = x[1]
			x2 = x[2]
			μ  = z[1]
			D  = u[1]
		#endregion

	#* Plant Equations	
    ##ODEs -
    dx[1] = x1  *(μ - D )
    dx[2] = D  *(x2_f - x2)          -   μ*x1/Y  

    ##Objective function - term inside integration
    dq = D*x1

    ##Alg Equations
    alg = μ - ( μ_max*x2/(km + x2 + k1*x2^2)   )

    return dx, dq, alg
end


#region->#generic function to Integrate a DAE using SunDials integrator
function Integrate_Plant_DAE(x0, q0, z0, u0, tspan)
        
		#region-> For Debugging (Arguments that are passed
			# x0 = [1.0; 1.0]
			# q0 = [0.0]
			# z0 = [0.53]
			# u0 = [0.1]
			# tspan = (0.0, 10.0)
		#endregion

		#region-> Loading Parameters
			Nx = size(x0,1)
			Nz = size(z0,1)
			Nu = size(u0,1)
			dx0  = 0*copy(x0)
			dq0  = 0*copy(q0)
			dz0  = 0*copy(z0)  #dummy variable analogous to xdot
			alg0 = 0*copy(z0)
		#endregion

	## function to collect variables and equations in correct order for passing to DAE integrator
	function fun(Uout,DU,U,p,t)

		x = U[1:Nx]    
		q = U[Nx+1]
		z = U[Nx+1+1 : Nx+1+Nz]
		u = p

		xdot = Plant_Eqns(x,z,u)[1]
		dq   = Plant_Eqns(x,z,u)[2]
		alg  = Plant_Eqns(x,z,u)[3]

		#Differential Equation for States (in form 0 = dx/dt - f(x,z,u))
		for nx = 1:Nx
		Uout[nx] = DU[nx] - xdot[nx]
		end

		#Differential Equation for Objective
		Uout[Nx+1] = DU[Nx+1] - dq

		#Algebraic Equations  (in form 0 = g(x,z,u))
		for nz = 1:Nz
		Uout[Nx+1+nz] = alg[nz]
		end

	end

	## Integrate the Plant for 1 time step (Solving DAE using Sundials)
		U₀  = vcat(x0,q0,z0)      	#Initial Guess for all the variables (Differential, Objective, Algebraic) in DAE
		DU₀ = vcat(dx0,dq0,alg0)  	#Initial Guess for their gradients (zero for the algebraic equations)
		id_differential_vars = vcat( ones(Nx+1), zeros(Nz))  #Identifier for which variables are differential and which are algebraic

		prob = DAEProblem(fun, DU₀,  U₀, tspan,  u0, differential_vars = id_differential_vars)

		#Solve the DAE
		sol = Sundials.solve(prob,IDA())
		
		xf = sol.u[end][1:Nx]
		qf = sol.u[end][Nx+1]
		zf = sol.u[end][Nx+1+1 : Nx+1+Nz]

	## Plot the Profile for all [x, q, z]
		plotlyjs()
		plot(sol)
  
	return xf, qf, zf
end
#endregion


## todo -> #generic function to integrate an ODE 
# function Integrate_Plant_ODE(x0, q0, u0, tspan)

# end

#Test