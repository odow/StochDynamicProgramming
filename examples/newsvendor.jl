#  Copyright 2015, Vincent Leclere, Francois Pacaud and Henri Gerard
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# Test SDDP with the newsvendor case study
#############################################################################

include("../src/SDDP.jl")
include("../src/simulate.jl")
include("../src/SDDPoptimize.jl")


using CPLEX
using JuMP

N_STAGES = 2
N_SCENARIOS = 1


function cost_t(x, u, w)
    h = .5
    p = 3

    cost = 0
    if x[1] >= 0
        cost += h*x[1]
    else
        cost += -p*x[1]
    end
    return cost + u[1]
end


function dynamic(x, u, w)
    return x + u -w
end


function init_problem()
    # Instantiate model:
    x0 = 4
    model = SDDP.LinearDynamicLinearCostSPmodel(N_STAGES, 1, 1, x0, cost_t, dynamic)
    solver = CplexSolver()
    params = SDDP.SDDPparameters(solver, N_SCENARIOS)

    return model, params
end


function solve_newsvendor()
    model, params = init_problem()
    V = optimize(model, params)
    law0 = NoiseLaw([0.],[1.])
    law1 = NoiseLaw([1., 2., 3., 4.], [.05, .75, .15, .05])
    law  = [law0,law1] 

    #aleas = simulate_scenarios(law ,(1, model.stageNumber, 1))
    aleas = simulate(law , 1)
    
    
    costs, stocks = forward_simulations(model, params, V, 1, aleas)
    println(stocks)
    println(costs)
    
    u0opt = law1.support[2]
    u1opt = law1.support - min(u0opt,law1.support)
    cost = u1opt*law1.proba + 0.5*u0opt + (2-5)*law1.support*law1.proba -2*min(u0opt,law1.support)*law1.proba
     println(u0opt)
     println(u1opt)
     println(cost)
    
end

@time solve_newsvendor()
