#=trying to solve and check Q4 on Problem Solving Coursework - 23/01/2022 =#

include("functionsmain.jl")
include("generalmain.jl")

regionNumber = 4

A_N = 0.6 #intial transfer coefficient A
B_N = 0.0 #intial transfer coefficient B

boundaries = zeros(Float64, regionNumber-1) #number of boundaries
potentialEnergy = zeros(Float64, regionNumber) #number of potentials

#=defining the variables in the question =#
a = 2.5053 *ang #the step potentials locations kinda thing
infty = 100*eV #approximating infitite potential energy as a large number(energy)
A = 0.6  #1.0 
q = 0.4 / ang #2/ang #0.0 * ang
Vcalc = ((hbar^2)/(2*mass)) * ((9*pi^2)/(16*a^2) + q^2) #expression for V0 I got...
n=3 #number of bound states
#Vcalc = ((hbar^2*pi^2)/(2*mass*(4*a^2))) * (n + 0.5)^2
V0 = 2*eV #5.0* eV #the height of the potential
println("Vcalc: ", Vcalc/eV, " eV")



#=setting up the boundaries=#
boundaries[1] = 0
boundaries[2] = a
boundaries[3] = 2*a
#boundaries[4] = 2*a

#=setting up the potential =#
potentialEnergy[1] = infty
potentialEnergy[2] = 0
potentialEnergy[3] = V0
potentialEnergy[4] = infty
#potentialEnergy[5] = infty






#=plot the potential well - region number representation=#
#potentialWellPlot = bar(1:regionNumber, potentialEnergy/eV, title = "Potential well (region number representation)", xlabel = "region number", ylabel = "E / eV", label = "", ylim = (0, V0/eV + 0.5))
#display(potentialWellPlot)

#=
#=defining the width of each bar=#
widthPlot = zeros(Float64, regionNumber)
widthPlot[1] = a
widthPlot[2] = boundaries[2] - boundaries[1]
widthPlot[3] = boundaries[3] - boundaries[2]
widthPlot[3] = a

#widthPlot[4] = boundaries[4] - boundaries[3]
#widthPlot[5] = a

#=new boundary values just for plotting and a better representation; plotting each bar in the middle of each region, nothing complicated =#
boundariesPlot = zeros(Float64, regionNumber)
boundariesPlot[1] = boundaries[1] - widthPlot[1]/2 #-2.5*a
boundariesPlot[2] = boundaries[2] - widthPlot[2]/2 #-1.5*a
boundariesPlot[3] = 0 
boundariesPlot[4] = boundaries[3] + widthPlot[4]/2 #1.5*a
boundariesPlot[5] = boundaries[4] + widthPlot[5]/2 #2.5*a



#=plot the potential well - region number representation=#
potentialWellPlot2 = bar(boundariesPlot/ang,potentialEnergy/eV, bar_width=widthPlot/ang, title = "Potential well (position representation)", xlabel = "x / ang", ylabel = "E / eV", label = "",  ylim = (0, V0/eV + 0.5))
display(potentialWellPlot2)

=#

#=range of energy values to calculate t11E over =#
energyStart = 0.01 * eV
energyEnd =  (V0+5) * eV #3.9800035035292347 * eV
energyStep = 0.01 * eV 

energyValues = energyStart:energyStep:energyEnd 

#=array to store t11(E)=#
t11E = zeros(ComplexF64, length(energyValues)) #t11(E)

#=calculate all t11=#
for e in 1:(length(energyValues))
    ABE2, kE2, T1NE2 = transferMatrix(mass, energyValues[e], potentialEnergy, boundaries, regionNumber, A_N, B_N)
    t11E[e] = T1NE2[1]
end

#=find all the roots of t11E=#
energyeigen = findRoots(real(t11E), energyValues/eV)

#=parameters of the graphs boundplot =#
xStart = boundaries[1] - 0.2*a
xEnd = boundaries[3] + 0.2*a
xStep = 0.01*ang


#=plot bound states =#
for e in 1:(length(energyeigen))
    caption = "Bound State: $e    Energy = $(energyeigen[e])eV"
    ABE2, kE2, T1NE2 = transferMatrix(mass, energyeigen[e]*eV, potentialEnergy, boundaries, regionNumber, A_N, B_N)
    boundplot, yvalues = plotWaveFunction(ABE2, kE2, regionNumber, xStep, boundaries, xStart, xEnd, caption)
    #boundplot = plot!(ylim = (-0.0000000000005, 0.0000000000005))
    display(boundplot)
    #savefig(boundplot, "boundplot$(replace(replace(string(Dates.now()), ":" => "_"), "."=> "_")).png")
end

#=the energy I calculated using the TISE =#
expectedEnergy = ((hbar^2) / (2*mass)) * ((9*pi^2)/(16*a^2))
println("Expected energy: ", expectedEnergy/eV, " eV")
V0calc2 = (expectedEnergy +  ((hbar^2) / (2*mass)) * q^2) / eV
println("V0calc2: ", V0calc2, " eV")
println("Calculated energy: ", energyeigen[2], " eV")


#=plotting the wave function in the question=#

    xvals1 =  (xStart:xStep:boundaries[1]) #/ ang
    xvals2 =  (boundaries[1]:xStep:boundaries[2]) #/ ang
    xvals3 =  (boundaries[2]:xStep:boundaries[3]) #/ ang
    xvals4 =  (boundaries[3]:xStep:boundaries[4]) #/ ang
    xvals5 =  (boundaries[4]:xStep:xEnd) #/ ang
    
    yvals1 = 0*xvals1
    yvals2 = -A*sinh.(q*(xvals2.+2*a)) # 2.3*sin.(2*(xvals2.+2*a)) 
    yvals3 = sin.( (3*pi.*xvals3) / (4*a) )
    yvals4 = A*sinh.(q*(2*a.-xvals4))
    yvals5 = 0*xvals5
    
    wavefunctionPlot = plot(title = "Wavefunction (from Q4)", xlabel = "x / ang", ylabel = "psi(x)", label = "") #, ylim = (0, V0/eV + 0.5)
    wavefunctionPlot = plot!(xvals1, yvals1,  linewidth = 3, label = "")  #, ylim = (0.001, 0.005) , label ="psi$j(x)"
    wavefunctionPlot = plot!(xvals2, yvals2,  linewidth = 3, label = "")
    wavefunctionPlot = plot!(xvals3, yvals3,  linewidth = 3, label = "")
    wavefunctionPlot = plot!(xvals4, yvals4,  linewidth = 3, label = "")
    wavefunctionPlot = plot!(xvals5, yvals5,  linewidth = 3, label = "")
    display(wavefunctionPlot)
    

#= 4(d) classically allowed probability =#
xvalsAll = [xvals1; xvals2; xvals3; xvals4; xvals5]
#yvalsAll = [yvals1; yvals2; yvals3; yvals4; yvals5]
#
#yvalsSquared = yvalsAll.^2

if length(energyeigen) >=2
    AB, k, T1N = transferMatrix(mass, energyeigen[2]*eV, potentialEnergy, boundaries, regionNumber, A_N, B_N)
    firstexcited, firstexcitedValues = plotWaveFunction(AB, k, regionNumber, xStep, boundaries, xStart, xEnd, "First energy eigen state")
end
yvalsSquared = firstexcitedValues.^2

wavefunctionSquaredPlot = plot(xvalsAll, real(yvalsSquared) ,title = "Wavefunction Squared ", xlabel = "x / ang", ylabel = "psi^2(x)", linewidth = 3, label = "") #, ylim = (0, V0/eV + 0.5)
display(wavefunctionSquaredPlot)