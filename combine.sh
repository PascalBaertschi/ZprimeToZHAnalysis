#!/bin/bash

# How to run:
#  source combine.sh -m alpha
# To filter jobs:
#  source combine.sh -m alpha XZhnnb_M

#option=""
#option="--freezeNuisanceGroups=theory --run=blind"
#for VBF categories alone:
#option="--freezeNuisanceGroups=theory --cminDefaultMinimizerStrategy 0"
#for combined categories:
#option="--rMin -3 --rMax 3 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_MaxCalls=999999999 --X-rtd MINIMIZER_analytic --cminDefaultMinimizerPrecision 1E-12"
#for all combined:
option="--rMin -8 --rMax 8 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_MaxCalls=999999999 --X-rtd MINIMIZER_analyti --cminDefaultMinimizerPrecision 1E-14"
#option="-H ProfileLikelihood"


if [[ "$1" != "-m" ]] || [[ "$2" == "" ]] || [[ "$3" == "" ]]
then
    echo Select a method and a filter:
    echo "  -m:          " $(ls datacards/)
    return
fi

higgsCombine() {
    card=$1
    method=$2
    inputfile=datacards/$method/$card.txt
    outputfile=combine/$method/$card.txt
    mass=$(echo $card | tr -dc '0-9')
    echo Running $method mass $mass on $inputfile...
    #echo $inputfile $outputfile $card
    > $outputfile
    if echo "$option" | grep -q "blind"; then
      echo '1' >> $outputfile
    fi

    #combine -M AsymptoticLimits --datacard $inputfile -m $mass $option
    combine -M AsymptoticLimits --datacard $inputfile -m $mass $option | grep -e Observed -e Expected | awk '{print $NF}' >> $outputfile
#    combine -M ProfileLikelihood --datacard $inputfile --significance --expectSignal=1 -m $mass $option | grep -e Significance -e value | awk '{print $NF}' | sed -e 's/)//g' >> $outputfile
#    combine -M MaxLikelihoodFit --datacard $inputfile -m $mass $option | grep Best | sed -e 's/\// /g' | awk '{print $4,$5,$6}' | sed -e 's/ /\n/g' >> $outputfile
#    combine -M MaxLikelihoodFit --datacard $inputfile --saveShapes --saveWithUncertainties -v 3 -n _$card | grep Best | sed -e 's/\// /g' | awk '{print $4,$5,$6}' | sed -e 's/ /\n/g' >> $outputfile # -t -1 generates errors # --robustFit=1 -m 0 --keepFailures
#    mv -f mlfit_$card.root combine/postfit/ # --plots
#    python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py --vtol=0.01 -f text combine/postfit/mlfit_$card.root | sed -e 's/!/ /g' -e 's/,/ /g' | tail -n +2 > combine/postfit/pull_$card.txt # -g rootfiles/pull_$1$2\_M$3.root
#    ##python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py --vtol=0.0000001 -f text  mlfit.root  | sed -e 's/!/ /g' -e 's/,/ /g' | tail -n +2 > pulls.txt
#    # Full CLs
#    outputfile=combine/$method/$card-fullCLs.txt
#    combine -M HybridNew --frequentist --testStat LHC -H ProfileLikelihood --datacard $inputfile -m $mass $option --fork 32 | grep -e 'Limit.*/' | awk '{print $4}' > $outputfile
#    combine -M HybridNew --frequentist --testStat LHC -H ProfileLikelihood --datacard $inputfile -m $mass --expectedFromGrid=0.025 $option --fork 32 | grep -e 'Limit.*/' | awk '{print $4}' >> $outputfile
#    combine -M HybridNew --frequentist --testStat LHC -H ProfileLikelihood --datacard $inputfile -m $mass --expectedFromGrid=0.16 $option --fork 32 | grep -e 'Limit.*/' | awk '{print $4}' >> $outputfile
#    combine -M HybridNew --frequentist --testStat LHC -H ProfileLikelihood --datacard $inputfile -m $mass --expectedFromGrid=0.5 $option --fork 32 | grep -e 'Limit.*/' | awk '{print $4}' >> $outputfile
#    combine -M HybridNew --frequentist --testStat LHC -H ProfileLikelihood --datacard $inputfile -m $mass --expectedFromGrid=0.84 $option --fork 32 | grep -e 'Limit.*/' | awk '{print $4}' >> $outputfile
#    combine -M HybridNew --frequentist --testStat LHC -H ProfileLikelihood --datacard $inputfile -m $mass --expectedFromGrid=0.975 $option --fork 32 | grep -e 'Limit.*/' | awk '{print $4}' >> $outputfile
}


for card in datacards/$2/$3*.txt
do
    analysis=$(basename $(dirname $card))
    signal=$(basename $card .txt)
    #echo Running with method $2 on $signal
    higgsCombine $signal $2 &
done

## Clean
wait
rm higgsCombine*.root
rm roostats-*
rm mlfit*.root

echo -e "\e[00;32mAll clear\e[00m"

# combine -M HybridNew --frequentist --datacard datacards/dijet/XVHah_M2600.txt --significance -H ProfileLikelihood --fork 32 -T 2000 -i 3

# python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/leeFromUpcrossings.py bands.root graph 1000 4500 --fit=best

# To derive Impacts (https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsWG/SWGuideNonStandardCombineUses#Nuisance_parameter_impacts):
# text2workspace.py datacards/dijet/XVHah_M2000.txt
# produces .root file in the same directory of the datacard
# combine -M MultiDimFit -n _initialFit_Test --algo singles --redefineSignalPOIs r --robustFit 1 -m 2000 -d XVHah_M2000.root
# combineTool.py -M Impacts -d XVHah_M2000.root -m 2000 --robustFit 1 --doFits --parallel 10
# combineTool.py -M Impacts -d XVHah_M2000.root -m 2000 -o impacts.json
# plotImpacts.py -i impacts.json -o impacts
