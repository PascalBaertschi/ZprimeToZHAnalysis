#!/bin/bash

#python signal.py -a -b

#python alpha.py -a -b -d | grep @ > alpha_norm.txt

# python alpha.py -c nnb -b -d &
# python alpha.py -c nnbb -b -d &
# wait


echo Waiting all background processes to finish...
wait

source combineCards.sh -m alpha

for cat in XZHnnb XZHnnbb XWHenb XWHenbb XWHmnb XWHmnbb XZHeeb XZHeebb XZHmmb XZHmmbb XZHnn XZHll XWHsl XZHsl XVHsl
do
    source combine.sh -m alpha "$cat"_M
    wait
    python limit.py -m alpha -c $cat -b
    echo Category $cat completed. Now waiting 60 seconds...
    sleep 60
done

#source combine.sh -m monoH monoHnn_M
#wait

#python limit.py -m alpha -c XZHnnb -b
#python limit.py -m alpha -c XZHnnbb -b
#python limit.py -m alpha -c XWHenb -b
#python limit.py -m alpha -c XWHenbb -b
#python limit.py -m alpha -c XWHmnb -b
#python limit.py -m alpha -c XWHmnbb -b
#python limit.py -m alpha -c XZHeeb -b
#python limit.py -m alpha -c XZHeebb -b
#python limit.py -m alpha -c XZHmmb -b
#python limit.py -m alpha -c XZHmmbb -b

#python limit.py -m alpha -c XZHnn -b
#python limit.py -m alpha -c XZHll -b

#python limit.py -m alpha -c XZHsl -b
#python limit.py -m alpha -c XWHsl -b
#python limit.py -m alpha -c XVHsl -b


echo -e "\e[00;32mAll clear\e[00m"


#source combineTest.sh datacards/alpha/XZHnnb_M1000.txt
#source combineTest.sh datacards/alpha/XZHnnb_M3000.txt
#source combineTest.sh datacards/alpha/XZHnnbb_M1000.txt
#source combineTest.sh datacards/alpha/XZHnnbb_M3000.txt
#source combineTest.sh datacards/alpha/XWHenb_M1000.txt
#source combineTest.sh datacards/alpha/XWHenb_M3000.txt
#source combineTest.sh datacards/alpha/XWHenbb_M1000.txt
#source combineTest.sh datacards/alpha/XWHenbb_M3000.txt
#source combineTest.sh datacards/alpha/XWHmnb_M1000.txt
#source combineTest.sh datacards/alpha/XWHmnb_M3000.txt
#source combineTest.sh datacards/alpha/XWHmnbb_M1000.txt
#source combineTest.sh datacards/alpha/XWHmnbb_M3000.txt
#source combineTest.sh datacards/alpha/XZHeeb_M1000.txt
#source combineTest.sh datacards/alpha/XZHeeb_M3000.txt
#source combineTest.sh datacards/alpha/XZHeebb_M1000.txt
#source combineTest.sh datacards/alpha/XZHeebb_M3000.txt
#source combineTest.sh datacards/alpha/XZHmmb_M1000.txt
#source combineTest.sh datacards/alpha/XZHmmb_M3000.txt
#source combineTest.sh datacards/alpha/XZHmmbb_M1000.txt
#source combineTest.sh datacards/alpha/XZHmmbb_M3000.txt

#source combineTest.sh datacards/alpha/XZHsl_M1000.txt
#source combineTest.sh datacards/alpha/XZHsl_M3000.txt
#source combineTest.sh datacards/alpha/XWHsl_M1000.txt
#source combineTest.sh datacards/alpha/XWHsl_M3000.txt
#source combineTest.sh datacards/alpha/XVHsl_M1000.txt
#source combineTest.sh datacards/alpha/XVHsl_M3000.txt

source combine.sh -m AZh AZh_M
wait
source combine.sh -m BBAZh BBAZh_M
wait
python limit.py -b -m 2HDM

echo -e "\e[00;32mAll clear\e[00m"
