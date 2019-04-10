#!/bin/bash 
function run {
	echo $2 ":"
	spaceex -g $1$2.cfg -m $1$3.xml -vl | grep 'Forbidden\|Computing reachable states done after'
}

function graphplot {
	echo $2 ":"
	spaceex -g $1$2.cfg -m $1$3.xml -o $1$2.gen -vl | grep 'Forbidden\|Computing reachable states done after'
	graph -Tpng --bitmap-size 1024x1024 -C -B -q0.5 $1$2.gen > $1$2.png
}

echo "ISS"
run SpaceStation/ ISSF01-ISS01 ISSF01
run SpaceStation/ ISSF01-ISU01 ISSF01
run SpaceStation/ ISSC01-ISS02 ISSC01
run SpaceStation/ ISSC01-ISU02 ISSC01
# Plot
echo "Producing plot for ISSF01-ISS02 over time horizon of 20"
graphplot SpaceStation/ ISSF01-ISS01-plot ISSF01-plot
echo "Producing plot for ISSC01-ISS02 over time horizon of 20"
graphplot SpaceStation/ ISSC01-ISS02-plot ISSC01-plot

echo "Spacecraft Rendezvous"
run Rendezvous/ SR01 SRNA01-SR0_
run Rendezvous/ SR02 SRA01-SR0_
# Plot
echo "Producing plot for SRNA01-SR0_ over time horizon of 120"
graphplot SpaceStation/ SR01-plot SRNA01-SR0_
echo "Producing plot for SRA01-SR0_ over time horizon of 120"
graphplot SpaceStation/ SR01-plot SRA01-SR0_.xml


echo "Building"
run Building/ BLDC01-BDS01 BLDC01
run Building/ BLDF01-BDS01 BLDF01
run Building/ BLDF01-BDU01 BLDF01
run Building/ BLDF01-BDU02 BLDF01
run Building/ BLDC01-BDU02 BLDC01
# Plot
echo "Producing plot for BLDF01-BDS01 over time horizon of 1"
graphplot Building/ BLDF01-BDS01-plot BLDF01 
echo "Producing plot for BLDF01-BDS01 over time horizon of 20"
graphplot Building/ BLDF01-BDS01-plot-20 BLDF01 

echo
echo "Platoon"
run Platoon/ PLAD01-BND42 PLAD01-BND
run Platoon/ PLAD01-BND30 PLAD01-BND
run Platoon/ PLAN01-UNB50 PLAN01-UNB

# Plot
echo "Producing plot for PLAD01-BND30 over time horizon of 20"
graphplot Platoon/ PLAD01-BND30-plot PLAD01-BND 

echo
echo "Gearbox"
run Gear/ GRBX01-MES01 GRBX01-MES01

# Plot
echo "Producing plot for GRBX01-MES01 over time horizon of 20"
graphplot Gear/ GRBX01-MES01-plot GRBX01-MES01


