include("Plant.jl")
include("OCP.jl")

##* Defining Simulation Parameters

T0_sim = 0
Tf_sim = 60
dt_sim = 1

Tf_ocp = 20
dt_ocp = 1

NFE_sim = convert(Int, (Tf_sim - T0_sim) / dt_sim)
NFE_ocp = convert(Int, (Tf_ocp) / dt_ocp)

yₛₚ = hcat(1.5032 * ones(1, convert(Int, NFE_sim / 3)),     0.9951 * ones(1,convert(Int, NFE_sim / 3)),    0 * ones(1,convert(Int, NFE_sim / 3)))
yₛₚ = hcat(yₛₚ,                                                                                             0 * ones(1,convert(Int, NFE_ocp)))

            #region-> storing parameters into named tuple
            Sim_par = (
                        T0_sim =  T0_sim, 
                        Tf_sim =  Tf_sim, 
                        dt_sim =  dt_sim, 
                        NFE_sim = NFE_sim,

                        Tf_ocp =  Tf_ocp,
                        dt_ocp =  dt_ocp,
                        NFE_ocp = NFE_ocp
                        )
            #endregion

# Starting Conditions for Simulation
uk₋₁    = [0.3]
xk₀     = [1.0; 1.0]
qk₀     = [0.0]
zk₀     = [0.0]

    #region ->Creating array for plotting
        x_plot = []
        q_plot = []
        z_plot = []
        u_plot = []

    # First Point for Plotting
        push!(x_plot, xk₀)
        push!(q_plot, qk₀)
        push!(z_plot, zk₀)
    #endregion


##* Simulate Plant with NMPC
for k in 1:NFE_sim
    # k = 21
    
    global uk₋₁, xk₀, qk₀, zk₀

    yₛₚk = yₛₚ[k:k - 1 + NFE_ocp]

    # Solve OCP
    uk₀ = Solve_OCP(xk₀, uk₋₁, yₛₚk)

    # Simulate Plant
    xk₁, qk₁, zk₁ = Integrate_Plant_DAE(xk₀, qk₀, zk₀, uk₀, (0.0, dt_sim))

    #region -> Updating values for next iteration and Plotting
    xk₀   = copy(xk₁)
    qk₀   = copy(qk₁)
    zk₀   = copy(zk₁)
    uk₋₁  = copy(uk₀)

    ##Results for Plotting
    push!(x_plot, copy(xk₀))
    push!(q_plot, copy(qk₀))
    push!(z_plot, copy(zk₀))
    push!(u_plot, copy(uk₀))

    #endregion
 
end

##* Plot Results
plotlyjs()
# gr()

t_plot = collect(T0_sim:dt_sim:Tf_sim)

p11 = plot( t_plot,[x_plot[nIter][1]  for nIter in 1:NFE_sim+1],            label="x1")
p11 = plot!(t_plot,[x_plot[nIter][2]  for nIter in 1:NFE_sim+1],            label="x2")
p11 = plot!(t_plot, yₛₚ[1:NFE_sim + 1],                                      label="y_sp",           linetype=:steppost, linestyle=:dash)

p12 = plot(t_plot[1:end - 1], [u_plot[nIter][1] for nIter in 1:NFE_sim],    label="u",              linetype=:steppost)

fig1 = plot(p11, p12, layout=(2, 1))



##
