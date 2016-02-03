#  Copyright 2014, Vincent Leclere, Francois Pacaud and Henri Gerard
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# Model and solve the One-Step One Alea problem in different settings
# - used to compute the optimal control (in forward phase / simulation)
# - used to compute the cuts in the Backward phase
#############################################################################

using JuMP
using SDDP

"""
Solve the Bellman equation at time t starting at state x under alea xi
with the current evaluation of Vt+1

The function solve
min_u current_cost(t,x,u,xi) + current_Bellman_Value_{t+1}(dynamic(t,x,u,xi))
and can return the optimal control and a subgradient of the value of the
problem with respect to the initial state x


Parameters:
- model (SPmodel)
    the stochastic problem we want to optimize

- param (SDDPparameters)
    the parameters of the SDDP algorithm

- V (bellmanFunctions)
    the current estimation of Bellman's functions

- t (int)
    time step at which the problem is solved

- xt (Array{Float})
    current starting state

- xi (Array{float})
    current noise value

- returnOptNextStage (Bool)
    return the optimal state at t+1

- returnOptcontrol (Bool)
    return the optimal control

- returnSubgradient (Bool)
    return the subgradient

- returnCost (Bool)
    return the value of the problem

TODO: update returns
TODO: add types in function parameters

"""
C = [0.5;2.0]

function solve_one_step_one_alea(model, #::SDDP.LinearDynamicLinearCostSPmodel,
                                 param, #::SDDP.SDDPparameters,
                                 V, #::Vector{SDDP.PolyhedralFunction},
                                 t, #::Int64,
                                 xt, #::Vector{Float64},
                                 xi) #::Vector{Float64},

    lambdas = V[t+1].lambdas
    betas = V[t+1].betas
    println(" ")
            println("lambdas")
            println(lambdas)   
            println(" ")
            
    println(" ")
            println("betas")
            println(betas)   
            println(" ")
    # TODO: factorize the definition of the model in PolyhedralFunction
    m = Model(solver=param.solver)
    @defVar(m, x)
    @defVar(m, 0 <= u <= 10)
    @defVar(m, alpha)
    
    if (t==2)
          @defVar(m, X)
    end

    @addConstraint(m, state_constraint, x .== xt)
    # TODO: implement cost function in PolyhedralFunction
    
    if (t==2)
         @addConstraint(m, X >= x)
         @addConstraint(m, X >= xi)
    end
    
    #@addConstraint(m, cost >= .5*x)
    #@addConstraint(m, cost >= -3*x)

    for i=1:V[t+1].numCuts
        @addConstraint(m, betas[i] + lambdas[i]*model.dynamics(x, u, xi) .<= alpha)
    end
     
    if (t==1)
          @setObjective(m, Min,C[1]*u  + alpha)
    end
    
    if (t==2)
          @setObjective(m, Min,(C[2]-5)*xi -C[2]*X + alpha)
    end
    
    println(" ")
            println("model")
            println(m)   
            println(" ")
    
    #@setObjective(m, Min, cost + alpha + u)

    status = solve(m)
    solved = (string(status) == "Optimal")

    if solved
        optimalControl = getValue(u)
        # Return object storing results:
        result = SDDP.NextStep(
                          model.dynamics(xt, optimalControl, xi),
                          [optimalControl],
                          getDual(state_constraint),
                          getObjectiveValue(m))
    else
        # If no solution is found, then return nothing
        result = nothing
    end

    return solved, result
end
