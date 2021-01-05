using Pkg
Pkg.update()
Pkg.add("StatsPlots")
Pkg.add("DataFrames")
Pkg.add("Colors")
Pkg.add("Showoff")
Pkg.add("Distributions")
Pkg.add("HypothesisTests")
Pkg.add("StatsPlots")
Pkg.add("Dates")
Pkg.add("CSV")
Pkg.add("Plots")
Pkg.add("PlotlyJS")

using StatsPlots
using Plots
using DataFrames
using Colors
using Showoff
using Distributions
using HypothesisTests
using CSV
using Dates

dataDir = "SimulationOutput"
N = 100 # Number of simulations
T = 120  # Number of periods per simulation

S = zeros(T,N) # 2 dim Array that will save S
U = zeros(T,N) # 2 dim Array that will save S
PlotFiles = ["A1", "A2", "A3", "A4",
                "B11", "B12", "B13", "B14",
                "B21", "B22", "B23", "B24",
                "B31", "B32", "B33", "B34",
                "B41", "B42", "B43", "B44",
                "C11", "C12", "C13", "C14",
                "C21", "C22", "C23", "C24",
                "C31", "C32", "C33", "C34",
                "C41", "C42", "C43", "C44",
                "D11", "D12", "D13", "D14",
                "D21", "D22", "D23", "D24",
                "D31", "D32", "D33", "D34",
                "D41", "D42", "D43", "D44",
                "E1", "E2", "E3", "E4",
                "F11", "F12", "F13", "F14", "F15",
                "F21", "F22", "F23", "F24", "F25",
                "F31", "F32", "F33", "F34", "F35",
                "F41", "F42", "F43", "F44", "F45",
                "G11", "G12", "G13", "G14",
                "G21", "G22", "G23", "G24",
                "G31", "G32", "G33", "G34",
                "G41", "G42", "G43", "G44",
                "H11", "H12", "H13", "H14",
                "H21", "H22", "H23", "H24",
                "H31", "H32", "H33", "H34",
                "H41", "H42", "H43", "H44",
                "I11", "I12", "I13", "I14",
                "I21", "I22", "I23", "I24",
                "I31", "I32", "I33", "I34",
                "I41", "I42", "I43", "I44",
                "J11", "J12", "J13", "J14", "J15",
                "J21", "J22", "J23", "J24", "J25",
                "J31", "J32", "J33", "J34", "J35",
                "J41", "J42", "J43", "J44", "J45"]

PlotFiles = ["ISU_A1", "ISU_A2", "ISU_A3", "ISU_A4", "ISU_B11",
            "ISU_B12", "ISU_B13", "ISU_B14", "ISU_B21",
            "ISU_B22", "ISU_B23", "ISU_B24", "ISU_B31",
            "ISU_B32", "ISU_B33", "ISU_B34"]

for file in PlotFiles
    # load S and U
    for i in 1:N
        df = DataFrame(CSV.File(joinpath(dataDir, string("$(i)_OUTPUT_Simulation_",file,"_120_Days.txt"))));
        for t in 1:T
            S[t,i] = df[t,:Susceptibles]
            U[t,i] = df[t,:Undetected]
        end
    end


    Ssummary = zeros(T,3) # S lower bound, mean, upper bound per period
    Usummary= zeros(T,3)  # U lower bound, mean, upper bound per period

    # Compute and save Ssummary and Usummary

    CI = 0.99 # confidence interval
    CIlb = (1.0-CI)/2
    CIub = 1.0 - (1.0-CI)/2
    for t in 1:T
        aux = quantile(S[t,:], [CIlb,CIub])
        Ssummary[t,1] = aux[1]          # S lower bound
        Ssummary[t,2] = mean(S[t,:])    # S mean
        Ssummary[t,3] = aux[2]          # S upper bound

        aux = quantile(U[t,:], [CIlb,CIub ])
        Usummary[t,1] = aux[1]          # U lower bound
        Usummary[t,2] = mean(U[t,:])    # U mean
        Usummary[t,3] = aux[2]          # U upper bound
    end


    titlefontsize = 14
    legendfontsize = 10
    tickfontsize = 10
    guidefontsize = 10
    linewidth = 3

    plotlyjs()
    # Create plot with CI as ribbons
    plot(Ssummary[:,2], label = "Susceptibles (S)"; ribbon = (Ssummary[:,2]-Ssummary[:,1], Ssummary[:,3]-Ssummary[:,2]), fillalpha = 0.2)
    plot!(Usummary[:,2], label = "New Infections (U)"; ribbon = (Usummary[:,2]-Usummary[:,1], Usummary[:,3]-Usummary[:,2]), fillalpha = 0.2)
    plot!(title = "Susceptibles (S) and New Infected (U)", xlabel = "Days Since Reopening", ylabel = "Number of Individuals")
    plot!(titlefontsize = titlefontsize, legendfontsize = legendfontsize, tickfontsize = tickfontsize, guidefontsize = tickfontsize)

    # Save the summary files in CSV
    fileString = joinpath("Plots", string("Susceptibles_Simulation_",file,"_120_Days.csv"))
    CSV.write(fileString, DataFrame(Ssummary))
    fileString = joinpath("Plots", string("Undetected_Simulation_",file,"_120_Days.csv"))
    CSV.write(fileString, DataFrame(Usummary))

    # https://docs.juliaplots.org/latest/generated/attributes_subplot/
    # Bottom left corner of legend is placed at (x,y). Symbol values: `:none`; `:best`; `:inline`; `:inside`; `:legend`; any valid combination of `:(outer ?)(top/bottom ?)(right/left ?)`, i.e.: `:top`, `:topright`, `:outerleft`, `:outerbottomright` ... (note: only some may be supported in each backend)
    plot!(legend = (0.75, 0.5))
    plot!(formatter = :plain)
    # Foreground color of the legend.
    # Colors names in https://github.com/JuliaGraphics/Colors.jl/blob/master/src/names_data.jl
    #plot!(foreground_color_legend=:lightgray)
    plot!(foreground_color_grid= :white)
    plot!(background_color_legend=:gray92)


    fileString = joinpath("Plots", string("plot_OUTPUT_Simulation_",file,"_120_Days.png"))
    savefig(fileString)

    fileString = joinpath("Plots", string("plot_OUTPUT_Simulation_",file,"_120_Days.pdf"))
    savefig(fileString)
end
