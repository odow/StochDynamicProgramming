#  Copyright 2015, Vincent Leclere, Francois Pacaud and Henri Gerard
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
#  the actual optimization function 
#
#############################################################################


include("../src/SDDP.jl")
println("SDDP.jl file included");

println(" ");
println("All files well included");
println(" ");

println("beginning of instanciation")
test = LinearDynamicLinearCostSPmodel(3,[2;2;2],[2;2;2],[1;1],costFunction,dynamic);

cut1 = PolyhedralFunction([1;1],[1 0;0 1]);
cut2 = PolyhedralFunction([2;2],[2 0;0 2]);
cut3 = PolyhedralFunction([3;3],[3 0;0 3]);

cut = PolyhedralFunction[cut1,cut2,cut3];

#stockTrajectories = Any[[1.0;1.0],[2.0;2.0],[3.0;3.0]];
stockTrajectories = ones(1,3,2);

stocks = zeros(1,3,2)
opt_control = zeros(1,3,2)

partest = SDDPparameters(GLPKSolverLP(),1,[0;20]);

alea = zeros(2,3);
println("end of instanciation")

println(" ");
println("Launch of function optimize");
optimize(test,partest,cut,alea);

println("                _    _")
println("               (_)  | |");
println(" _____  _   _   _   | |");
println("|  _  || | | | | |  | |");
println("| |_| || |_| | | |  \\ /");
println("|_____||_____| |_|   o ");
