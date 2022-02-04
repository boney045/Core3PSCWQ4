
using Distributions
using Plots
using LinearAlgebra


function transferMatrix(m, E, U, x, regionNum, A, B)

    Coefficient = zeros(ComplexF64, regionNum, 2) #a matrix of size[regionNum x 2] to store transfer coefficients (As and Bs)
    Tmatrix = zeros(ComplexF64, regionNum-1, 4) #a matrix of size[regionNum x 4] to store transfer matrix coefficients
    T1N = [1.0+0*im 0.0+0*im; 0.0+0.0*im 1.0+0.0*im] #Transmission coefficient 1,N
    k = zeros(ComplexF64, regionNum) #an array to store wave factors
    Coefficient[regionNum, 1] = A #set the transfer coefficients in the matrix
    Coefficient[regionNum, 2] = B #set the transfer coefficients in the matrix
    Tmat = [0.0 0.0; 0.0 0.0]
    for i in (regionNum-1):-1:1
        k[i] = sqrt(Complex(2*m*(E-U[i])))/hbar
        k[i+1] = sqrt(Complex(2*m*(E-U[i+1])))/hbar

        tempMat = [Coefficient[i+1,1];  Coefficient[i+1,2]]

        if k[i] == 0
            Tmat = [(1-1im*k[i+1]*x[i])*exp(1im*k[i+1]*x[i]) (1+1im*k[i+1]*x[i])*exp(-1im*k[i+1]*x[i]); 1im*k[i+1]*exp(1im*k[i+1]*x[i]) -1im*k[i+1]*exp(-1im*k[i+1]*x[i])]
        elseif k[i+1] == 0
            Tmat = [(1/2)*exp(-1im*k[i]*x[i]) (1/(2*k[i]))*(k[i]*x[i]-1im)*exp(-1im*k[i]*x[i]); (1/2)*exp(1im*k[i]*x[i]) (1/(2*k[i]))*(k[i]*x[i]+1im)*exp(1im*k[i]*x[i])]
        else
            Tmat = 1/(2*k[i]) * [ (k[i]+k[i+1])*exp(-1im*x[i]*(k[i]-k[i+1])) (k[i]-k[i+1])*exp(-1im*x[i]*(k[i]+k[i+1])) ;
            (k[i]-k[i+1])*exp(1im*x[i]*(k[i]+k[i+1])) (k[i]+k[i+1])*exp(1im*x[i]*(k[i]-k[i+1]))] 
        end

        Tmatrix[i,1] = Tmat[1,1]
        Tmatrix[i,2] = Tmat[1,2]
        Tmatrix[i,3] = Tmat[2,1]
        Tmatrix[i,4] = Tmat[2,2]

        tempTransfer = Tmat * tempMat

        Coefficient[i, 1] = tempTransfer[1]
        Coefficient[i, 2] = tempTransfer[2]
    end

    for i in 1:1:(regionNum-1)
        t = [Tmatrix[i,1] Tmatrix[i,2]; Tmatrix[i,3] Tmatrix[i, 4]]
        T1N = T1N * t 
    end
    #println("Coefficient $Coefficient")
    return Coefficient, k, T1N
    
end


function plotWaveFunction(AB, k, regionNum, dx, boundary, startPos, endPos, caption)
    graph = plot(xlabel = "x / m", ylabel = "real(PSI)", title = caption)
    xvalues = startPos:dx:boundary[1]
    yvalues = []
    for j in 1:regionNum
        if j>1
            if j < regionNum
                xvalues = boundary[j-1]:dx:boundary[j]
            else
                xvalues = boundary[j-1]:dx:endPos
            end
        end
        psi = AB[j,1]*exp.(1im*k[j].*xvalues) + AB[j,2]*exp.(-1im*k[j].*xvalues)
        #println(psi)
        yvalues = vcat(yvalues, psi)
        graph = plot!(xvalues, real(psi),  linewidth = 3, label = "")  #, ylim = (0.001, 0.005) , label ="psi$j(x)"
    end
    return graph, yvalues
end

function setBoundary(boundary, startX, deltaX, regionNum)
    boundary[1] = startX
    dx = deltaX/(regionNum-2)
    for i in 1:regionNum-2
        boundary[i+1] = boundary[i] + dx
    end
    return boundary
end

function probTransmissionAndReflection(transmissionProb, reflectionProb, mass, energyValues, potentialEnergy, boundaries, regionNumber, A_N, B_N)
    #=calculate transmision and reflectoin probablity for each energy value =#
    for e in 1:(length(energyValues))
        #=call the transferMat function to get coefficients, k-values and T_1N for each energy in energyValues=#
        coefficients, k, T_1N = transferMatrix(mass, energyValues[e], potentialEnergy, boundaries, regionNumber, A_N, B_N )
        #=calculate Transmission and relfection for each energy and store it in TrE and RE=#
        transmissionProb[e] = abs(abs(1/T_1N[1,1])^2 * (k[regionNumber]/k[1])) 
        reflectionProb[e] = abs(T_1N[2,1]/T_1N[1,1])^2
    end
    return transmissionProb, reflectionProb
end


function equalTriangleBarrier(U, h, initalU, regionNum)
    U[1] = initalU
    U[regionNum] = initalU
    du = h/round(((regionNum-2)/2)+0.1)
    rem = regionNum%2
    for i in 2:regionNum-1
        if i <= round(regionNum/2 + 0.1) 
            U[i] = U[i-1] + du
        elseif rem==0 && i==regionNum/2 + 1
            U[i] = U[i-1] 
        else 
            U[i] = U[i-1] - du
        end
    end
    #U[regionNum] = initalU
    return U
end


function gaussiaBarrier(U, h, w, boundary, deltaX, regionNum)
    dx = deltaX/(regionNum-2) #boundary width 
    b= boundary[1] + deltaX/2 + dx/2 #boundary[Int(round(regionNum/2))]  # peak center

    c = w / (2*sqrt(2*log(2)))
    #println("c: ", c)

    for i in 1:regionNum-1
        U[i] = h*exp(-(boundary[i]-b)^2/(2*c^2))
        #println("U[$i]: ", U[i])
    end
    U[regionNum] = h*exp(-(boundary[regionNum-1]+dx-b)^2/(2*c^2))

    return U    
end

function harmonicPotential(U, k, boundary, regionNum, deltaX)
    dx = deltaX/(regionNum-2) #boundary width 
    shiftright = deltaX/2
    for i in 1:regionNum-1
        U[i] = (k*(boundary[i] - shiftright)^2)/2
    end
    U[regionNum] = (k*(boundary[regionNumber-1]+dx)^2)/2

    return U      
end

#see if you can make it general for n number of wells right now I am going to hard code - 24/12/2021
function potentialWell(U, heightU, boundary, startB, deltaX, regionNum)
    if regionNum != 3
        println("Region number is ", regionNum, ".")
        println("Region number is required to be 3 for a simple potential well")
        return -1
    end

    U[1] = heightU
    U[2] = 0
    U[3] = heightU

    boundary[1] = startB
    boundary[2] = boundary[1] + deltaX

    return U, boundary
end

function potentialWell2(U, heightU, boundary, startB, width, bwidth, regionNum)

    for i in 1:regionNum
        if mod(i,2) == 0
            U[i] = 0
        else
            U[i] = heightU
        end
    end

    boundary[1] = startB
    for b in 2:regionNum-1
        if mod(b,2) == 0
            boundary[b] = boundary[b-1] + width
        else
            boundary[b] = boundary[b-1] + bwidth
        end
    end
    return U, boundary
end


function findRoots(yvalues, xvalues)
    roots = Vector{Float64}()
    for i in 1:length(yvalues)-1
        current = yvalues[i]
        next = yvalues[i+1]
        if current == 0
            #println("root: ",current,  ", at index: ", i)
            #println("xvalues: ", xvalues[i] )
            push!(roots, round((xvalues[i]), sigdigits = 5))
            println()
        elseif current  * next < 0 
            #println("root between ", current, " and ", next, " at index: ", i," and ",i+1)
            #println("between xvalues: ", xvalues[i], " and ", xvalues[i+1] )
            push!(roots, round((xvalues[i]+xvalues[i+1])/2, sigdigits = 5))
        end
    end
    if yvalues[length(yvalues)] == 0
        #println("root: ",yvalues[length(yvalues)],  ", at index: ", length(yvalues))
        #println("xvalues: ", xvalues[length(yvalues)] )
        push!(roots, round(xvalues[length(yvalues)], sigdigits = 5))
    end
    return roots
end


#testarr = [-4,-3,-2,-1,1,2,3,4,3,4,3,2,1,0,-1,1,0]

#findRoots(testarr)