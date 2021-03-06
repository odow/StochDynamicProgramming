#  Copyright 2014, Vincent Leclere, Francois Pacaud and Henri Gerard
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# Define all types used in this module.
#############################################################################

abstract SPModel

type LinearDynamicLinearCostSPmodel <: SPModel
    # problem dimension
    stageNumber::Int64
    dimControls::Int64
    dimStates::Int64
    dimNoises::Int64

    # Bounds of states and controls:
    xlim
    ulim

    initialState

    costFunctions
    dynamics
    noises

    # TODO: add this attributes to model
    # lowerbounds#::Tuple{Vector{Float64}}
    # upperbounds#::Tuple{Vector{Float64}}
end



type SDDPparameters
    solver
    forwardPassNumber::Int64 # number of simulated scenario in the forward pass

    # TODO: add this attributes to SDDPparameters
    # initialization #TODO
    # stoppingTest #TODO
    # maxItNumber #TODO
end

type PolyhedralFunction
    #function defined by max_k betas[k] + lambdas[k,:]*x
    betas::Vector{Float64}
    lambdas::Array{Float64,2} #lambdas[k,:] is the subgradient

    # number of cuts:
    numCuts::Int64
end


type NextStep
    next_state::Array{Float64, 1}
    optimal_control::Array{Float64, 1}
    sub_gradient
    cost
end
