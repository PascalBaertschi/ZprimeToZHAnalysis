#!/bin/bash

if [[ "$1" != "-m" ]] || [[ "$2" == "" ]]
then
    echo Select a method:
    echo "  -m:          " $(ls datacards/)
    return
fi

cd datacards/$2/


# Step 1: merge 1 and 2 b-tag categories
if [[ "$2" == "alpha" ]]
then
#    for channel in XZHnn XWHen XWHmn XZHee XZHmm
#    do
#        echo Merging b-tag categories in channel $channel...
#        for card in "$channel"bb_M*.txt
#        do
#            mass="${card//[!0-9]/}"
#            combineCards.py "$channel"b="$channel"b_M$mass.txt "$channel"bb="$channel"bb_M$mass.txt > "$channel"_M$mass.txt
#        done
#    done

    # Step 2: merge lepton categories
    echo Merging 0b and bb categories...
    for card in XZHmmbb_M*.txt
    do
        mass="${card//[!0-9]/}"
	combineCards.py nn0b=XZHnn0b_M$mass.txt nnbb=XZHnnbb_M$mass.txt > XZHnn_M$mass.txt
	combineCards.py ee0b=XZHee0b_M$mass.txt eebb=XZHeebb_M$mass.txt > XZHee_M$mass.txt
	combineCards.py mm0b=XZHmm0b_M$mass.txt mmbb=XZHmmbb_M$mass.txt > XZHmm_M$mass.txt
	combineCards.py ee0b=XZHee0b_M$mass.txt eebb=XZHeebb_M$mass.txt mm0b=XZHmm0b_M$mass.txt mmbb=XZHmmbb_M$mass.txt > XZHll_M$mass.txt
	combineCards.py nn0b=XZHnn0b_M$mass.txt ee0b=XZHee0b_M$mass.txt mm0b=XZHmm0b_M$mass.txt > XZHsl0b_M$mass.txt
	combineCards.py nnbb=XZHnnbb_M$mass.txt eebb=XZHeebb_M$mass.txt mmbb=XZHmmbb_M$mass.txt > XZHslbb_M$mass.txt
	combineCards.py nn0b=XZHnn0b_M$mass.txt ee0b=XZHee0b_M$mass.txt mm0b=XZHmm0b_M$mass.txt  nnbb=XZHnnbb_M$mass.txt eebb=XZHeebb_M$mass.txt mmbb=XZHmmbb_M$mass.txt > XZHsl_M$mass.txt
	combineCards.py nn0bVBF=XZHVBFnn0bVBF_M$mass.txt nnbbVBF=XZHVBFnnbbVBF_M$mass.txt > XZHnnVBF_M$mass.txt
	combineCards.py ee0bVBF=XZHVBFee0bVBF_M$mass.txt eebbVBF=XZHVBFeebbVBF_M$mass.txt > XZHeeVBF_M$mass.txt
	combineCards.py mm0bVBF=XZHVBFmm0bVBF_M$mass.txt mmbbVBF=XZHVBFmmbbVBF_M$mass.txt > XZHmmVBF_M$mass.txt
	combineCards.py ee0bVBF=XZHVBFee0bVBF_M$mass.txt eebbVBF=XZHVBFeebbVBF_M$mass.txt mm0bVBF=XZHVBFmm0bVBF_M$mass.txt mmbbVBF=XZHVBFmmbbVBF_M$mass.txt > XZHllVBF_M$mass.txt
	combineCards.py nn0bVBF=XZHVBFnn0bVBF_M$mass.txt ee0bVBF=XZHVBFee0bVBF_M$mass.txt mm0bVBF=XZHVBFmm0bVBF_M$mass.txt > XZHsl0bVBF_M$mass.txt
	combineCards.py nnbbVBF=XZHVBFnnbbVBF_M$mass.txt eebbVBF=XZHVBFeebbVBF_M$mass.txt mmbbVBF=XZHVBFmmbbVBF_M$mass.txt > XZHslbbVBF_M$mass.txt
	combineCards.py nn0bVBF=XZHVBFnn0bVBF_M$mass.txt ee0bVBF=XZHVBFee0bVBF_M$mass.txt mm0bVBF=XZHVBFmm0bVBF_M$mass.txt  nnbbVBF=XZHVBFnnbbVBF_M$mass.txt eebbVBF=XZHVBFeebbVBF_M$mass.txt mmbbVBF=XZHVBFmmbbVBF_M$mass.txt > XZHslVBF_M$mass.txt
	combineCards.py nn0b=XZHnn0b_M$mass.txt nn0bVBF=XZHVBFnn0bVBF_M$mass.txt nnbb=XZHnnbb_M$mass.txt nnbbVBF=XZHVBFnnbbVBF_M$mass.txt> XZHnncomb_M$mass.txt
	combineCards.py ee0b=XZHee0b_M$mass.txt ee0bVBF=XZHVBFee0bVBF_M$mass.txt eebb=XZHeebb_M$mass.txt eebbVBF=XZHVBFeebbVBF_M$mass.txt> XZHeecomb_M$mass.txt
	combineCards.py mm0b=XZHmm0b_M$mass.txt mm0bVBF=XZHVBFmm0bVBF_M$mass.txt mmbb=XZHmmbb_M$mass.txt mmbbVBF=XZHVBFmmbbVBF_M$mass.txt> XZHmmcomb_M$mass.txt
	combineCards.py ee0b=XZHee0b_M$mass.txt ee0bVBF=XZHVBFee0bVBF_M$mass.txt eebb=XZHeebb_M$mass.txt eebbVBF=XZHVBFeebbVBF_M$mass.txt mm0b=XZHmm0b_M$mass.txt mm0bVBF=XZHVBFmm0bVBF_M$mass.txt mmbb=XZHmmbb_M$mass.txt mmbbVBF=XZHVBFmmbbVBF_M$mass.txt> XZHllcomb_M$mass.txt
	combineCards.py nn0b=XZHnn0b_M$mass.txt nn0bVBF=XZHVBFnn0bVBF_M$mass.txt ee0b=XZHee0b_M$mass.txt ee0bVBF=XZHVBFee0bVBF_M$mass.txt mm0b=XZHmm0b_M$mass.txt mm0bVBF=XZHVBFmm0bVBF_M$mass.txt> XZHsl0bcomb_M$mass.txt
	combineCards.py nnbb=XZHnnbb_M$mass.txt nnbbVBF=XZHVBFnnbbVBF_M$mass.txt eebb=XZHeebb_M$mass.txt eebbVBF=XZHVBFeebbVBF_M$mass.txt mmbb=XZHmmbb_M$mass.txt mmbbVBF=XZHVBFmmbbVBF_M$mass.txt> XZHslbbcomb_M$mass.txt
	combineCards.py nn0b=XZHnn0b_M$mass.txt nn0bVBF=XZHVBFnn0bVBF_M$mass.txt ee0b=XZHee0b_M$mass.txt ee0bVBF=XZHVBFee0bVBF_M$mass.txt mm0b=XZHmm0b_M$mass.txt  mm0bVBF=XZHVBFmm0bVBF_M$mass.txt nnbb=XZHnnbb_M$mass.txt nnbbVBF=XZHVBFnnbbVBF_M$mass.txt eebb=XZHeebb_M$mass.txt eebbVBF=XZHVBFeebbVBF_M$mass.txt mmbb=XZHmmbb_M$mass.txt mmbbVBF=XZHVBFmmbbVBF_M$mass.txt> XZHslcomb_M$mass.txt
    done
fi


# Step 3: merge all hadronic
if [[ "$2" == "dijet" ]]
then
    echo Merging all hadronic categories...
    echo -n "Signal "
    for signal in XWH XZH XVH
        do
        echo -n " $signal"
        for card in XWHwrhpb_M*.txt
        do
            mass="${card//[!0-9]/}"
            combineCards.py wrhpb="$signal"wrhpb_M$mass.txt wrhpbb="$signal"wrhpbb_M$mass.txt wrlpb="$signal"wrlpb_M$mass.txt wrlpbb="$signal"wrlpbb_M$mass.txt zrhpb="$signal"zrhpb_M$mass.txt zrhpbb="$signal"zrhpbb_M$mass.txt zrlpb="$signal"zrlpb_M$mass.txt zrlpbb="$signal"zrlpbb_M$mass.txt > "$signal"ah_M$mass.txt
            #combineCards.py wrhpb="$signal"wrhpb_M$mass.txt wrhpbb="$signal"wrhpbb_M$mass.txt wrlpb="$signal"wrlpb_M$mass.txt wrlpbb="$signal"wrlpbb_M$mass.txt zrhpb="$signal"zrhpb_M$mass.txt zrhpbb="$signal"zrhpbb_M$mass.txt zrlpb="$signal"zrlpb_M$mass.txt zrlpbb="$signal"zrlpbb_M$mass.txt wrhphp="$signal"wrhphp_M$mass.txt wrhplp="$signal"wrhplp_M$mass.txt wrlphp="$signal"wrlphp_M$mass.txt zrhphp="$signal"zrhphp_M$mass.txt zrhplp="$signal"zrhplp_M$mass.txt zrlphp="$signal"zrlphp_M$mass.txt > "$signal"ah_M$mass.txt #wrlplp="$signal"wrlplp_M$mass.txt zrlplp="$signal"zrlplp_M$mass.txt 
        done
    done
    echo
fi


if [[ "$2" == "AZh" ]] || [[ "$2" == "BBAZh" ]]
then
    echo Copying cards...
    for mA in {800..2000..100}
    do
        for cat in nnb nnbb eeb eebb mmb mmbb
        do
            old=XZH"$cat"_M"$mA"
            new="$2$cat"_M"$mA"
            new2="$2$cat"_M"$((mA-50))"
            # copy card
            cp ../alpha/$old.txt $new.txt
            cp ../alpha/$old.txt $new2.txt
            sed -i "s/$old/$new/g" $new.txt
            sed -i "s/XZH$cat/$2$cat/g" $new.txt
            sed -i "s/$old/$new2/g" $new2.txt
            sed -i "s/XZH$cat/$2$cat/g" $new2.txt
        done
    done
    echo Combining cards...
    for mA in {800..2000..50}
    do
        # combine cards
        combineCards.py nnb="$2"nnb_M"$mA".txt nnbb="$2"nnbb_M"$mA".txt > "$2"nn_M"$mA".txt
        combineCards.py eeb="$2"eeb_M"$mA".txt eebb="$2"eebb_M"$mA".txt mmb="$2"mmb_M"$mA".txt mmbb="$2"mmbb_M"$mA".txt > "$2"ll_M"$mA".txt
        combineCards.py nnb="$2"nnb_M"$mA".txt nnbb="$2"nnbb_M"$mA".txt eeb="$2"eeb_M"$mA".txt eebb="$2"eebb_M"$mA".txt mmb="$2"mmb_M"$mA".txt mmbb="$2"mmbb_M"$mA".txt > "$2"_M"$mA".txt
    done
fi


if [[ "$2" == "monoH" ]]
then
    echo Converting cards...
    for mZ in {800..4000..50}
    do
        for mA in {300..1200..10}
        do
            echo "Merging cards for $mZ $mA" 
            for cat in nnb nnbb
            do
                oldcard=XZH"$cat"_M"$mZ".txt
                newcard=monoH"$cat"_MZ"$mZ"_MA"$mA".txt
                oldwork=XZH"$cat".root
                newwork=monoH"$cat".root
                oldsignal=XZH"$cat"_M"$mZ"
                newsignal=monoH"$cat"_MZ"$mZ"_MA"$mA"
#                # check if workspace exists
#                if [ ! -f "../../workspace/$newwork" ]
#                then
#                    echo Mass point mZ $mZ and mA $mA does not exist
#                    continue
#                fi
                # copy card
                if [ -f "../alpha/$oldcard" ]
                then
                    cp ../alpha/$oldcard $newcard
                    sed -i "s/"$oldsignal"/"$newsignal"/g" $newcard
                    sed -i "s/"$oldwork"/"$newwork"/g" $newcard
#                    sed -i "s/b_eig0                 param     0.0                 1.0/b_eig0                 param     0.0                 1.0/g" $newcard
#                    sed -i "s/b_eig1                 param     0.0                 1.0/b_eig1                 param     0.0                 1.0/g" $newcard
#                    sed -i "s/b_eig2                 param     0.0                 1.0/b_eig2                 param     0.0                 1.0/g" $newcard
#                    sed -i "s/b_eig3                 param     0.0                 1.0/b_eig3                 param     0.0                 1.0/g" $newcard
#                    sed -i "s/b_eig4                 param     0.0                 1.0/b_eig4                 param     0.0                 1.0/g" $newcard
#                    sed -i "s/b_eig5                 param     0.0                 1.0/b_eig5                 param     0.0                 1.0/g" $newcard
                else
                    oldcard2=XZH"$cat"_M"$((mZ-50))".txt
                    oldsignal2=XZH"$cat"_M"$((mZ-50))"
                    cp ../alpha/$oldcard2 $newcard
                    sed -i "s/"$oldsignal2"/"$newsignal"/g" $newcard
                    sed -i "s/"$oldwork"/"$newwork"/g" $newcard
                fi
            done
            # combine cards
            if [ -f "monoHnnb_MZ"$mZ"_MA"$mA".txt" ] && [ -f "monoHnnbb_MZ"$mZ"_MA"$mA".txt" ]
            then
                combineCards.py nnb=monoHnnb_MZ"$mZ"_MA"$mA".txt nnbb=monoHnnbb_MZ"$mZ"_MA"$mA".txt > monoHnn_MZ"$mZ"_MA"$mA".txt
            fi
        done
    done
fi

if [[ "$2" == "combo" ]]
then
    echo Copying datacards...
    cp ../alpha/* .
    cp ../dijet/* .
    echo Merging Alpha and Dijet...
    for signal in XWH XZH XVH
        do
        echo -n " $signal"
        for card in XVHsl_M*.txt
        do
            mass="${card//[!0-9]/}"
            # Below 1 TeV there is no ah
            if [ ! -f "$signal"ah_M$mass.txt ]
            then
                combineCards.py nnb="$signal"nnb_M$mass.txt nnbb="$signal"nnbb_M$mass.txt enb="$signal"enb_M$mass.txt enbb="$signal"enbb_M$mass.txt mnb="$signal"mnb_M$mass.txt mnbb="$signal"mnbb_M$mass.txt eeb="$signal"eeb_M$mass.txt eebb="$signal"eebb_M$mass.txt mmb="$signal"mmb_M$mass.txt mmbb="$signal"mmbb_M$mass.txt > "$signal"_M$mass.txt
#                combineCards.py nnb="$signal"nnb_M$mass.txt nnbb="$signal"nnbb_M$mass.txt enb="$signal"enb_M$mass.txt enbb="$signal"enbb_M$mass.txt mnb="$signal"mnb_M$mass.txt mnbb="$signal"mnbb_M$mass.txt eeb="$signal"eeb_M$mass.txt eebb="$signal"eebb_M$mass.txt mmb="$signal"mmb_M$mass.txt mmbb="$signal"mmbb_M$mass.txt > "$signal"_M$mass.txt
            else
                combineCards.py nnb="$signal"nnb_M$mass.txt nnbb="$signal"nnbb_M$mass.txt enb="$signal"enb_M$mass.txt enbb="$signal"enbb_M$mass.txt mnb="$signal"mnb_M$mass.txt mnbb="$signal"mnbb_M$mass.txt eeb="$signal"eeb_M$mass.txt eebb="$signal"eebb_M$mass.txt mmb="$signal"mmb_M$mass.txt mmbb="$signal"mmbb_M$mass.txt wrhpb="$signal"wrhpb_M$mass.txt wrhpbb="$signal"wrhpbb_M$mass.txt wrlpb="$signal"wrlpb_M$mass.txt wrlpbb="$signal"wrlpbb_M$mass.txt zrhpb="$signal"zrhpb_M$mass.txt zrhpbb="$signal"zrhpbb_M$mass.txt zrlpb="$signal"zrlpb_M$mass.txt zrlpbb="$signal"zrlpbb_M$mass.txt > "$signal"_M$mass.txt
#                combineCards.py wrhpb="$signal"wrhpb_M$mass.txt nnbb="$signal"nnbb_M$mass.txt > "$signal"_M$mass.txt
            fi
        done
    done
    echo
fi
cd ../..

echo -e "\e[00;32mAll clear\e[00m"
