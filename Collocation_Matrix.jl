


function collocation_matrix(NCP, Poly)
    
    # NCP = 3
    # Poly = "Radau"

    #todo - error handling properly implemented. 
    # struct incomplete_code <: Exception end



    #*Radau Polynomials
    if Poly == "Radau"

        if NCP == 2
            # t[k,i] = [0, 1]

            Pdot = [ -1.0  1.0]

            P = [0  1]

        elseif NCP == 3
            #t[k,i] = [0, 0.333333, 1]

            Pdot = [   -2.0   1.5   0.499999		
                        2.0  -4.5   2.5		]

            P = [-4.44E-16  0   1] 

        elseif NCP == 4
            #t[k,1] = [0, 0.155051, 0.644949, 1]

            Pdot = [    -4.13939   3.22475   1.16784   -0.253197		
                        1.73939  -3.56784   0.775255   1.0532		
                        -3.0       5.53197  -7.53197    5.0]	

            P = [-1.78E-15  0   8.88E-16 1]

        else
            # throw(incomplete_code())
        end



    #*Legendre Polynomials
    elseif Poly == "Legendre"
        #todo - Complete matrices for Legendre
        if NCP == 2
            # throw(incomplete_code())
        elseif NCP == 3
            # throw(incomplete_code())
        elseif NCP == 4
            # throw(incomplete_code())
        else
            # throw(incomplete_code())
        end
    end






    return Pdot, P 
end



# myPdot, myP = collocation_matrix(4, "Radau")
# myPdot
# myP