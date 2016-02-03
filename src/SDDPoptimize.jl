#  Copyright 2015, Vincent Leclere, Francois Pacaud and Henri Gerard
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
#  the actual optimization function
#
#############################################################################


include("forwardBackwardIterations.jl")
include("utility.jl")
include("simulate.jl")

"""
TODO: add docstring
TODO: move initialize in proper module
TODO: fix initialize

"""
function get_null_value_functions_array(model::SDDP.SPModel)

    V = Vector{SDDP.PolyhedralFunction}(model.stageNumber)
    for t = 1:model.stageNumber
        V[t] = get_null_value_functions()
    end

    return V
end


function initialize_value_functions( model::SDDP.LinearDynamicLinearCostSPmodel,
                                     param::SDDP.SDDPparameters,
                                     law#::NoiseLaw
                        )

    V_null = get_null_value_functions_array(model)
    V = Array{SDDP.PolyhedralFunction}(model.stageNumber)

    #aleas = simulate_scenarios(law,
    #                          (param.forwardPassNumber,
    #                            model.stageNumber, 1))
    aleas = simulate(law, 1)

    n = param.forwardPassNumber

    V[end] = SDDP.PolyhedralFunction([-7.4;-10.25;-10.85;-11], [-4. 0.; -2. 0.; 0. 0.;1. 0.], 4)

    #stockTrajectories = forward_simulations(model,
     #                   param,
      #                  V_null,
       #                 n,
        #                aleas)[2]

    stockTrajectories = zeros(1,2,1)
    stockTrajectories[1,1,1] = 4.0
    stockTrajectories[1,1,1] = 0.0

    backward_pass(model,
                  param,
                  V,
                  stockTrajectories,
                  law,
                  true)
    return V
end



"""
Make a forward pass of the algorithm

Simulate a scenario of noise and compute an optimal trajectory on this
scenario according to the current value functions.

Parameters:
- model (SPmodel)
    the stochastic problem we want to optimize

- param (SDDPparameters)
    the parameters of the SDDP algorithm


Returns :
- V::Array{PolyhedralFunction}
    the collection of approximation of the bellman functions

"""
function optimize(model::SDDP.SPModel,
                  param::SDDP.SDDPparameters)

    # Initialize value functions:
    #law = NoiseLaw([0., 1., 2., 3.], [.2, .4, .3, .1])
    law0 = NoiseLaw([0.],[1.])
    law1 = NoiseLaw([0., 1., 2., 3.], [.25, .25, .25, .25])
    law  = [law0,law1]
    V = initialize_value_functions(model, param, law)
    #aleas = simulate_scenarios(law ,(param.forwardPassNumber, model.stageNumber, 1))
    aleas = simulate(law, 1)
    stopping_test::Bool = false
    iteration_count::Int64 = 0

    n = param.forwardPassNumber

    for i = 1:2
        stockTrajectories = forward_simulations(model,
                            param,
                            V,
                            n,
                            aleas)[2]
        backward_pass(model,
                      param,
                      V,
                      stockTrajectories,
                      law)
        # TODO: stopping test
          
          println(stockTrajectories)
        iteration_count+=1;
    end

    return V
end
