#!/bin/bash

cd tests/params
mkdir scan
for ealpha in -2 -1.5 -1; do
	for egamma in -2 -1.5 -1; do
		for ewd0 in -2 -1.5 -1 -0.5; do
			fname=ea$ealpha"_"eg$egamma"_"ewd$ewd0
			echo $fname
			sed    's/ eWD_alpha .*/ eWD_alpha     '$ealpha'/' p_base.ini > p_$fname.ini
			sed -i 's/ eWD_gamma .*/ eWD_gamma     '$egamma'/' p_$fname.ini
			sed -i 's/ eWD .*/ eWD     '$ewd0'/' p_$fname.ini
			sed -i 's/exptName .*/exptName        par_'$fname'/' p_$fname.ini
			mv p_$fname.ini scan/p_$fname.ini
		done
	done
done

cd ../..
mkdir logs
for ealpha in -2 -1.5 -1; do
	for egamma in -2 -1.5 -1; do
		for ewd0 in -2 -1.5 -1 -0.5; do
			fname=ea$ealpha"_"eg$egamma"_"ewd$ewd0
			echo $fname
			./tests/pf.test tests/params/scan/p_$fname.ini > logs/out_$fname.log &
		done
	done
done

