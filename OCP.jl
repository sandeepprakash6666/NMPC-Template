#=
28/8/20 - Exercise 5 from MPC modeule 2019 version at NTNU.
Creates collocation matrix of 3rd order Radau polynomial
Solves OCP taking in initiial differential state, previous MV (for penalizing input moves) and setpoint trajectory
=#
using Plots
using JuMP
using Ipopt

include("parameters.jl")

function Collocation_Matrix()
      #Radau
      t1 = 0.155051
      t2 = 0.644949
      t3 = 1.0

      M1 = [
            t1 1 / 2 * t1^2 1 / 3 * t1^3
            t2 1 / 2 * t2^2 1 / 3 * t2^3
            t3 1 / 2 * t3^2 1 / 3 * t3^3
            ]
      M2 = [
            1 t1 t1^2
            1 t2 t2^2
            1 t3 t3^2
            ]

      M = M1 * inv(M2)
      return M
end

##*
function Solve_OCP(x0, u₋₁, y_sp)

            #region -> Value of Arguments for Debugging
                  # x0 = [1.0; 1.0]
                  # u₋₁ = 0.3
                  # y_sp = append!(1.5302*ones(10), 0.9951*ones(10) )
            #endregion

            #region-> Loading parameters
                  #Model parameters
                  model_par = Model_par()
                        x2_f  = model_par[:x2_f]
                        km    = model_par[:km]
                        k1    = model_par[:k1]
                        Y     = model_par[:Y]
                        μ_max = model_par[:μ_max]
                        
                  #generic code (variable Parameters)
                  Bounds = Var_bounds()
                        Nx = Bounds[:Nx]
                        Nz = Bounds[:Nz]
                        Nu = Bounds[:Nu]
                        x_guess = Bounds[:x_guess]
                        z_guess = Bounds[:z_guess]
                        u_guess = Bounds[:u_guess]
                        us_x = Bounds[:us_x]
                        ls_x = Bounds[:ls_x]
                        us_z = Bounds[:us_z]
                        ls_z = Bounds[:ls_z]
                        us_u = Bounds[:us_u]
                        ls_u = Bounds[:ls_u]
            #endregion


      ##*OCP Parameters
      T0  =                                                     0.0
      Tf  = haskey(Sim_par,   :Tf_ocp)    ? Sim_par.Tf_ocp    : 30
      dt  = haskey(Sim_par,   :dt_ocp)    ? Sim_par.dt_ocp    : 1
      NFE = haskey(Sim_par,   :NFE_ocp)   ? Sim_par.NFE_ocp   : 30
      NCP = 3

      Solve_OCP     = true
      Display_Plots = false

                  #region -> Setting Initial guesses for quadrature and dx
                  q0   = 0
                  dq0  = 0
                  dx0  = 0*copy(x_guess)
                  alg0 = 0*copy(z_guess) 
                  #endregion

      ##* Defining Solver
      model1 = Model(with_optimizer(Ipopt.Optimizer))

            #region-> #*Define Variables and Objective
                  ## Declare Variables - Scaled between 0-1
                  @variable(model1, x[1:Nx, 1:NFE, 1:NCP])
                  @variable(model1, z[1:Nz, 1:NFE, 1:NCP])
                  @variable(model1, u[1:Nu, 1:NFE])

                  #unscaled variables
                  @variable(model1, q[    1,       1:NFE, 1:NCP])
                  @variable(model1, dx_us[1:Nx,    1:NFE, 1:NCP])
                  @variable(model1, dq[   1,       1:NFE, 1:NCP])
                  @variable(model1, alg[  1:Nz,    1:NFE, 1:NCP])


                  #region -> Set Variable Bounds AND Initial Guesses 
                  for nx in 1:Nx, nz in 1:Nz, nu in 1:Nu, nfe in 1:NFE, ncp in 1:NCP     
                        set_lower_bound(x[nx, nfe, ncp], 0)
                        # set_upper_bound(x[1 , nfe, ncp], 1)

                        set_lower_bound(z[nz, nfe, ncp], 0)
                        #set_upper_bound(z[nz, nfe, ncp], 999)

                        #set_lower_bound(dx[nx, nfe, ncp], 0)
                        #set_upper_bound(dx[nx, nfe, ncp], 999)

                        #set_lower_bound(alg_coll[nz, nfe, ncp], 0)
                        #set_upper_bound(alg_coll[nz, nfe, ncp], 999)

                        set_lower_bound(u[nu, nfe], 0)
                        set_upper_bound(u[nu, nfe], 1)
                  end

                  #Set Initial Guesses
                  for nx in 1:Nx, nz in 1:Nz, nu in 1:Nu, nfe in 1:NFE, ncp in 1:NCP
                        set_start_value(x[nx, nfe, ncp],    x_guess[nx])
                        set_start_value(z[nz, nfe, ncp],    z_guess[nz])
                        set_start_value(dx_us[nx, nfe, ncp],dx0[nx])
                        set_start_value(alg[nz, nfe, ncp],  alg0[nz])
                        set_start_value(u[nu, nfe],         u_guess[nu])
                        set_start_value(q[1, nfe, ncp],     q0)
                        set_start_value(dq[1, nfe, ncp],    dq0)
                  end
                  #endregion

                  #region-> Expressions for Unscaling Variables (makes it easier to write DAE Equation) #?add variables and change indices 
                  @NLexpressions(model1, begin

                        x1[nfe in 1:NFE, ncp in 1:NCP],           x[1, nfe, ncp]            *  (us_x[1]     - ls_x[1])   + ls_x[1]
                        x2[nfe in 1:NFE, ncp in 1:NCP],           x[2, nfe, ncp]            *  (us_x[2]     - ls_x[2])   + ls_x[2]

                        μ[nfe in 1:NFE, ncp in 1:NCP],            z[1, nfe, ncp]            *  (us_z[1]     - ls_z[1])   + ls_z[1]

                        D[nfe in 1:NFE],                          u[1, nfe]                 *  (us_u[1]     - ls_u[1])   + ls_u[1] 

                  end)
                  #endregion

            ## Objective
            @NLobjective(model1, Min,  sum( (x1[nfe,NCP] - y_sp[nfe])^2   for nfe in 1:NFE )  + 0.05*(D[1] - u₋₁[1])^2 + sum( 0.05*(D[nfe] - D[nfe-1])^2 for nfe in 2:NFE)              )
            
            #endregion

            #region-> #*Define Constraints

            @NLconstraints(model1, begin
                  #Defining the model ODEs in each line
                  Constr_ODE1[nfe in 1:NFE, ncp in 1:NCP], dx_us[1, nfe, ncp]      == x1[nfe, ncp] * (μ[nfe, ncp] - D[nfe])
                  Constr_ODE2[nfe in 1:NFE, ncp in 1:NCP], dx_us[2, nfe, ncp]      == D[nfe]       * (x2_f - x2[nfe, ncp])           - μ[nfe, ncp] * x1[nfe, ncp] / Y
                  #In case of more states - pattern
                  #Constr_ODE999[nfe=1:NFE, ncp=1:NCP], dx[999,nfe,ncp] ==
            end)

            @NLconstraints(model1, begin
                  #Collocation for Quadrature/ Objective Function
                  #t = 0
                  Constr_quad_dot0[nfe = 1    , ncp in 1:NCP], dq[1, nfe, ncp]  == (x1[nfe,ncp] - y_sp[nfe])^2 + 0.5*(D[nfe]  - u₋₁[1]    )^2
                  #t = 1..N
                  Constr_quad_dot[nfe in 2:NFE, ncp in 1:NCP], dq[1, nfe, ncp]  == (x1[nfe,ncp] - y_sp[nfe])^2 + 0.5*(D[nfe] - D[nfe-1])^2
                  #In case of objective depends on end state only -> dq = dx[2] etc
            end)


            @NLconstraints(model1, begin
                  #Defining Model Algebraic Equations in each line
                  Constr_Alg1[nfe in 1:NFE, ncp in 1:NCP], μ[nfe, ncp] == (  μ_max * x2[nfe, ncp] / (km + x2[nfe, ncp] + k1 * x2[nfe, ncp]^2) )
                  #In case of more states - pattern
                  #Constr_Alg999[nfe=1:NFE, ncp=1:NCP], alg[999,nfe,ncp] ==
            end)

            @NLconstraints(model1, begin
                  #Defining Inequality Constraints in each line
                  Constr_Ineq1[nfe in 1:1   ], D[nfe] - u₋₁[1]       <= 0.05
                  Constr_Ineq2[nfe in 2:NFE ], D[nfe] - D[nfe-1]     <= 0.05
                  #In case of more states - pattern
                  #Constr_Ineq999[nfe=1:NFE, ncp=1:NCP], alg[999,nfe,ncp] ==
            end)


                  #region-> #generic code -> Collocation Equation for Differential Equations AND Objective Function
                        ## Creating a Radau collocation Matrix for NCP = 3
                        collMat = Collocation_Matrix()
                        @NLconstraints(model1, begin
                              #Collocation Equation 
                              #t = 0
                              Constr_Coll_Diff0[nx in 1:Nx,  nfe = 1, ncp in 1:NCP],        x[nx, nfe, ncp]    == (x0[nx] -ls_x[nx])/(us_x[nx] - ls_x[nx])              + dt * sum(collMat[ncp, i] * dx_us[nx, nfe, i] for i = 1:NCP)/(us_x[nx] - ls_x[nx])
                              Constr_Coll_quad0[             nfe = 1, ncp in 1:NCP],        q[1, nfe, ncp]     == q0                                                    + dt * sum(collMat[ncp, i] * dq[1, nfe, i]  for i = 1:NCP)
                              #t = 1 ... (N-1)
                              Constr_Coll_Diff[nx in 1:Nx,   nfe in 2:NFE, ncp = 1:NCP],    x[nx, nfe, ncp]     == x[nx, nfe-1, NCP]                                    + dt * sum(collMat[ncp, i] * dx_us[nx, nfe, i] for i = 1:NCP)/(us_x[nx] - ls_x[nx])
                              Constr_Coll_quad[              nfe in 2:NFE, ncp = 1:NCP],    q[1, nfe, ncp]      == q[1, nfe-1, NCP]                                     + dt * sum(collMat[ncp, i] * dq[1, nfe, i]  for i = 1:NCP)
                        end)
                  #endregion
            #endregion

      ##*Solve the model
      if Solve_OCP == true

            optimize!(model1)
            JuMP.termination_status(model1)
            JuMP.solve_time(model1::Model)

            #get results into Julia variables
                  #scaled values
                  star_x = JuMP.value.(x[:, :, NCP])
                  star_x = cat(x0, star_x, dims = 2)
                  star_u = JuMP.value.(u)
                  star_z = JuMP.value.(z[:, :, NCP])
                  #Unscaled variables
                  star_x1 = JuMP.value.(x1[:, NCP])
                  star_x2 = JuMP.value.(x2[:, NCP])
                        star_x1 = cat(x0[1], star_x1, dims = 1)     
                        star_x2 = cat(x0[2], star_x2, dims = 1)     
                  star_D = JuMP.value.(D)

      end

            #region-> #*Plotting Solution
            if Display_Plots == true

                  t_plot = collect(T0:dt:Tf)    #Returns NFE+1 dimensional vector
                  #
                  #choose backend for plots
                  plotlyjs()
                  #
                  p11 = plot(t_plot, star_x1,               label = "x1")
                  p11 = plot!(t_plot, star_x[2, :],         label = "x2")
                  p11 = plot!(t_plot[1:end-1], y_sp[:],     label = "y_sp",   linetype = :steppost, linestyle = :dash)
                  p12 = plot(t_plot[1:end-1], star_D,       label = "D",      linetype = :steppost)
                  
                  fig1 = plot(p11, p12, layout = (2, 1))
            end
            #endregion

      star_MPC = (star_u[:,1] .* (us_u - ls_u)) .+ ls_u

      return star_MPC
end

