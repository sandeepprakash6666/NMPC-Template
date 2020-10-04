
function Var_bounds()
    
    ##? Unscaled Initial POints/ Guesses
    x0_us = [1.5;  1.0]
    z0_us = [0.53]
    u0_us = [0.3]

    ##? Lower and Upper scales
    ls_x = [0.0;   0.0]
    us_x = [10.0;  10.0]

    ls_z = [0.0]
    us_z = [1.0]

    ls_u = [0.0]
    us_u = [1.0]


    #region-> Collecting to Dictionary
        myDict = Dict( 
        
        :ls_x => ls_x,
        :us_x => us_x, 

        :ls_z => ls_z, 
        :us_z => us_z, 

        :ls_u => ls_u, 
        :us_u => us_u,

        #Dimensions
        :Nx => size(x0_us, 1),
        :Nz => size(z0_us, 1),
        :Nu => size(u0_us, 1),

        # ##Scaled Points
        :x_guess =>    (x0_us - ls_x)      ./ (us_x - ls_x),
        :z_guess =>    (z0_us - ls_z)      ./ (us_z - ls_z),
        :u_guess =>    (u0_us - ls_u)      ./ (us_u - ls_u)
        )
    #endregion

   return myDict
end




function Model_par()
    
    myDict = Dict(

        :x2_f   => 4,
        :km     => 0.12,
        :k1     => 0.4545,
        :Y      => 0.4,
        :μ_max  => 0.53
    )

    return myDict
end


function Plant_par()
    
    myDict = Dict(

        :x2_f   => 4,
        :km     => 0.12,
        :k1     => 0.4545,
        :Y      => 0.4,
        :μ_max  => 0.53
    )

    return myDict
end


