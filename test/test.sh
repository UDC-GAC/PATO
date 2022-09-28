#!/bin/bash

export OMP_NUM_THREADS=4
FLAVOUR=gnu

# Set up
rm -rf output
mkdir output

[ $# -gt 0 ] && FLAVOUR=$1

PATO=../target/$FLAVOUR/pato.serial
PATO_PAR=../target/$FLAVOUR/pato.release

TFO_FILE=input/tfo.fa
TTS_FILE=input/tts.fa

FAILED=0

echo -e "\033[0;31m!! Warning:\033[0m tests run using the GNU flavour of the tool by default"

# Test 0
printf "Test #0 -- Default run... "
$PATO $TFO_FILE $TTS_FILE output/foo0 1> /dev/null
$PATO_PAR $TFO_FILE $TTS_FILE output/bar0 1> /dev/null

sort output/foo0.out > output/fuu0.out
sort output/bar0.out > output/ber0.out
sort output/foo0.summary > output/fuu0.summary
sort output/bar0.summary > output/ber0.summary
sort ref/0.ref.out > output/0.out
sort ref/0.ref.summary > output/0.summary

diff output/fuu0.out output/0.out > output/0.diff
diff output/ber0.out output/0.out > output/1.diff
diff output/fuu0.summary output/0.summary > output/2.diff
diff output/ber0.summary output/0.summary > output/3.diff

RESULT=0
for i in 0 1 2 3; do
    NUM_DIFFS=$(wc -l < output/$i.diff)
    RESULT=$(($RESULT + $NUM_DIFFS))
done

if [[ $RESULT -eq 0 ]]
then
    printf "\033[0;32mOK\033[0m\n"
else
    printf "\033[0;31mFAILED\033[0m\n"
    FAILED=$((FAILED + 1))
fi

# Test 1
printf "Test #1 -- Default run + Triplex extended format... "
$PATO -of 1 $TFO_FILE $TTS_FILE output/foo1 1> /dev/null
$PATO_PAR -of 1 $TFO_FILE $TTS_FILE output/bar1 1> /dev/null

sort output/foo1.out > output/fuu1.out
sort output/bar1.out > output/ber1.out
sort output/foo1.summary > output/fuu1.summary
sort output/bar1.summary > output/ber1.summary
sort ref/1.ref.out > output/1.out
sort ref/1.ref.summary > output/1.summary

diff output/fuu1.out output/1.out > output/0.diff
diff output/ber1.out output/1.out > output/1.diff
diff output/fuu1.summary output/1.summary > output/2.diff
diff output/ber1.summary output/1.summary > output/3.diff

RESULT=0
for i in 0 1 2 3; do
    NUM_DIFFS=$(wc -l < output/$i.diff)
    RESULT=$(($RESULT + $NUM_DIFFS))
done

if [[ $RESULT -eq 0 ]]
then
    printf "\033[0;32mOK\033[0m\n"
else
    printf "\033[0;31mFAILED\033[0m\n"
    FAILED=$((FAILED + 1))
fi

# Test 2
printf "Test #2 -- Default run + Purine strand error reference format... "
$PATO -er 1 $TFO_FILE $TTS_FILE output/foo2 1> /dev/null
$PATO_PAR -er 1 $TFO_FILE $TTS_FILE output/bar2 1> /dev/null

sort output/foo2.out > output/fuu2.out
sort output/bar2.out > output/ber2.out
sort output/foo2.summary > output/fuu2.summary
sort output/bar2.summary > output/ber2.summary
sort ref/2.ref.out > output/2.out
sort ref/2.ref.summary > output/2.summary

diff output/fuu2.out output/2.out > output/0.diff
diff output/ber2.out output/2.out > output/1.diff
diff output/fuu2.summary output/2.summary > output/2.diff
diff output/ber2.summary output/2.summary > output/3.diff

RESULT=0
for i in 0 1 2 3; do
    NUM_DIFFS=$(wc -l < output/$i.diff)
    RESULT=$(($RESULT + $NUM_DIFFS))
done

if [[ $RESULT -eq 0 ]]
then
    printf "\033[0;32mOK\033[0m\n"
else
    printf "\033[0;31mFAILED\033[0m\n"
    FAILED=$((FAILED + 1))
fi

# Test 3
printf "Test #3 -- Default run + Third strand error reference format... "
$PATO -er 2 $TFO_FILE $TTS_FILE output/foo3 1> /dev/null
$PATO_PAR -er 2 $TFO_FILE $TTS_FILE output/bar3 1> /dev/null

sort output/foo3.out > output/fuu3.out
sort output/bar3.out > output/ber3.out
sort output/foo3.summary > output/fuu3.summary
sort output/bar3.summary > output/ber3.summary
sort ref/3.ref.out > output/3.out
sort ref/3.ref.summary > output/3.summary

diff output/fuu3.out output/3.out > output/0.diff
diff output/ber3.out output/3.out > output/1.diff
diff output/fuu3.summary output/3.summary > output/2.diff
diff output/ber3.summary output/3.summary > output/3.diff

RESULT=0
for i in 0 1 2 3; do
    NUM_DIFFS=$(wc -l < output/$i.diff)
    RESULT=$(($RESULT + $NUM_DIFFS))
done

if [[ $RESULT -eq 0 ]]
then
    printf "\033[0;32mOK\033[0m\n"
else
    printf "\033[0;31mFAILED\033[0m\n"
    FAILED=$((FAILED + 1))
fi

# Clean up
rm -rf output

# Check if failed
if [[ $FAILED -gt 0 ]]
then
    exit 1
fi
exit 0
