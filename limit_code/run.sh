

bkgtype=("air")
sigtype=("const")
obstype=("obs")
DMmass=("0p01" "0p05" "0p1" "0p5" "1" "5" "10" "50" "100")

mypath=$PWD
echo $mypath

for (( i=0; i<${#sigtype[@]}; i++ ))
do
  for (( j=0; j<${#bkgtype[@]}; j++ ))
  do
    for (( m=0; m<${#obstype[@]}; m++ ))
    do
      echo "Signal velocity: ${sigtype[$i]}, background material: ${bkgtype[$j]}, observation type: ${obstype[$m]}"
      if [ -d "UL_${sigtype[$i]}_${bkgtype[$j]}_${obstype[$m]}" ]; then
        rm -rf UL_${sigtype[$i]}_${bkgtype[$j]}_${obstype[$m]}
      fi
      mkdir UL_${sigtype[$i]}_${bkgtype[$j]}_${obstype[$m]}
      cd UL_${sigtype[$i]}_${bkgtype[$j]}_${obstype[$m]}
      for (( k=0; k<${#DMmass[@]}; k++ ))
      do
        echo "Dark Matter mass: ${DMmass[$k]}"
        cp -r ../datacard.txt datacard_${sigtype[$i]}_${DMmass[$k]}_${bkgtype[$j]}_${obstype[$m]}.txt
        sed -i "s/SIG/${sigtype[$i]}_${DMmass[$k]}/g" datacard_${sigtype[$i]}_${DMmass[$k]}_${bkgtype[$j]}_${obstype[$m]}.txt
        sed -i "s/BKG/${bkgtype[$j]}/g" datacard_${sigtype[$i]}_${DMmass[$k]}_${bkgtype[$j]}_${obstype[$m]}.txt
        sed -i "s/OBS/${obstype[$m]}/g" datacard_${sigtype[$i]}_${DMmass[$k]}_${bkgtype[$j]}_${obstype[$m]}.txt
        #(
        #text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO verbose --PO 'map=.*/h_sig:r_sig[1, 0, 10]' --PO 'map=.*/h_bkg:r_bkg[1, 0, 10]' datacard_${sigtype[$i]}_${DMmass[$k]}_${bkgtype[$j]}_${obstype[$m]}.txt
        #combine -M MultiDimFit -d datacard_${sigtype[$i]}_${DMmass[$k]}_${bkgtype[$j]}_${obstype[$m]}.root --setParameters r_sig=1,r_bkg=1 --X-rtd FAST_VERTICAL_MORPH --algo grid --points=50 --floatOtherPOIs=1 --redefineSignalPOIs r_sig,r_bkg -P r_sig --rMin 0 --rMax 10 -n ${sigtype[$i]}_${DMmass[$k]}_${bkgtype[$j]}_${obstype[$m]}_scan_r_sig
        #combine -M MultiDimFit -d datacard_${sigtype[$i]}_${DMmass[$k]}_${bkgtype[$j]}_${obstype[$m]}.root --setParameters r_sig=1,r_bkg=1 --X-rtd FAST_VERTICAL_MORPH --algo grid --points=50 --floatOtherPOIs=1 --redefineSignalPOIs r_sig,r_bkg -P r_bkg --rMin 0 --rMax 10 -n ${sigtype[$i]}_${DMmass[$k]}_${bkgtype[$j]}_${obstype[$m]}_scan_r_bkg
        #combine -M AsymptoticLimits datacard_${sigtype[$i]}_${DMmass[$k]}_${bkgtype[$j]}_${obstype[$m]}.root -n ${sigtype[$i]}_${DMmass[$k]}_${bkgtype[$j]}_${obstype[$m]} --redefineSignalPOIs r_sig,r_bkg --setParameters r_sig=1,r_bkg=1
        #combine -M AsymptoticLimits datacard_${sigtype[$i]}_${DMmass[$k]}_${bkgtype[$j]}_${obstype[$m]}.txt -n ${sigtype[$i]}_${DMmass[$k]}_${bkgtype[$j]}_${obstype[$m]}
        #) 2>&1 | tee log_${sigtype[$i]}_${DMmass[$k]}_${bkgtype[$j]}_${obstype[$m]} &
        #break
	combine -M AsymptoticLimits datacard_${sigtype[$i]}_${DMmass[$k]}_${bkgtype[$j]}_${obstype[$m]}.txt -n ${sigtype[$i]}_${DMmass[$k]}_${bkgtype[$j]}_${obstype[$m]} > log_${sigtype[$i]}_${DMmass[$k]}_${bkgtype[$j]}_${obstype[$m]} &
      done
      cd $mypath
    done
  done
done
#wait
