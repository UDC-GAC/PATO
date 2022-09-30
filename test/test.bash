#!/bin/bash

# Set up
exec 2> /dev/null

dir=$(dirname -- "$( readlink -f -- "$0"; )";)

flavour=gnu
[ $# -gt 0 ] && flavour=$1

rm -rf $dir/output
mkdir $dir/output

pato=$dir/../target/$flavour/pato.serial
pato_par=$dir/../target/$flavour/pato.release

tfo_file=$dir/input/tfo.fa
tts_file=$dir/input/tts.fa

failed=0
total=0

printf "\e[1mPATO: PArallel TriplexatOr\e[0m -- Test suite\n\n"

# TFO tests
declare -a tfo_opts=("" "-l 12 -L 34" "-l 12 -L 34 -e 10" "-l 14 -L 28 -e 9 -E 5" "-c 3 -g 5 -G 50" "-b 3" "-a true" "-e 2 -mrl 4" "-e 2 -mrl 7 -mrp 1" "-e 10 -dd 1 -dc 2" "-e 10 -dd 2 -dc 2" "-e 10 -dd 1 -dc 2 -dl true" "-e 10 -dd 2 -dc 2 -dl true")
declare -a tfo_flags=("" "-po true " "-mf true " "-fr false " "-mf true -fr false " "-ssd off ")

fails=0
counter=0
printf "         \e[1mTFO search\e[0m\n"
for opts in "${tfo_opts[@]}"; do
    for flags in "${tfo_flags[@]}"; do
        printf "\e[1m[  --  ]\e[0m Test #${counter} ${flags}${opts}"

        sort $dir/ref/tfo${counter}.out > $dir/output/tfo${counter}.ref.sorted.out
        sort $dir/ref/tfo${counter}.summary > $dir/output/tfo${counter}.ref.sorted.summary

        $pato $opts $flags -ss $tfo_file -o $dir/output/tfo${counter}.seq 1> /dev/null
        sort $dir/output/tfo${counter}.seq.out > $dir/output/tfo${counter}.seq.sorted.out
        sort $dir/output/tfo${counter}.seq.summary > $dir/output/tfo${counter}.seq.sorted.summary

        $pato_par $opts $flags -ss $tfo_file -o $dir/output/tfo${counter}.par 1> /dev/null
        sort $dir/output/tfo${counter}.par.out > $dir/output/tfo${counter}.par.sorted.out
        sort $dir/output/tfo${counter}.par.summary > $dir/output/tfo${counter}.par.sorted.summary

        diff $dir/output/tfo${counter}.seq.sorted.out $dir/output/tfo${counter}.ref.sorted.out > $dir/output/tfo${counter}.seq.diff.out
        diff $dir/output/tfo${counter}.par.sorted.out $dir/output/tfo${counter}.ref.sorted.out > $dir/output/tfo${counter}.par.diff.out
        diff $dir/output/tfo${counter}.seq.sorted.summary $dir/output/tfo${counter}.ref.sorted.summary > $dir/output/tfo${counter}.seq.diff.summary
        diff $dir/output/tfo${counter}.par.sorted.summary $dir/output/tfo${counter}.ref.sorted.summary > $dir/output/tfo${counter}.par.diff.summary

        result=0
        result=$(($result + $(wc -l < $dir/output/tfo${counter}.seq.diff.out)))
        result=$(($result + $(wc -l < $dir/output/tfo${counter}.par.diff.out)))
        result=$(($result + $(wc -l < $dir/output/tfo${counter}.seq.diff.summary)))
        result=$(($result + $(wc -l < $dir/output/tfo${counter}.par.diff.summary)))

        if [[ $result -eq 0 ]]; then
            printf "\r\e[1m[  \033[0;32mOK\033[0m  \e[1m]\e[0m\n"
        else
            printf "\r\e[1m[ \033[0;31mFAIL\033[0m \e[1m]\e[0m\n"
            failed=$(($failed + 1))
            fails=$((fails + 1))
        fi

        counter=$(($counter + 1))
    done
done
total=$(($total + $counter))
printf "         \e[1mRun:\e[0m ${counter} \e[1mOK\e[0m: $(($counter - $fails)) \e[1mFAIL:\e[0m ${fails}\n\n"

# TTS tests
declare -a tts_opts=("" "-l 12 -L 34" "-l 12 -L 34 -e 10" "-l 14 -L 28 -e 9 -E 5" "-c 3 -g 5 -G 50" "-b 3" "-a true" "-e 2 -mrl 4" "-e 2 -mrl 7 -mrp 1" "-e 10 -dd 1 -dc 2" "-e 10 -dd 2 -dc 2" "-e 10 -dd 1 -dc 2 -dl true" "-e 10 -dd 2 -dc 2 -dl true")
declare -a tts_flags=("" "-po true " "-mf true " "-fr off " "-mf true -fr off " "-ssd off ")

fails=0
counter=0
printf "         \e[1mTTS search\e[0m\n"
for opts in "${tts_opts[@]}"; do
    for flags in "${tts_flags[@]}"; do
        printf "\e[1m[  --  ]\e[0m Test #${counter} ${flags}${opts}"

        sort $dir/ref/tts${counter}.out > $dir/output/tts${counter}.ref.sorted.out
        sort $dir/ref/tts${counter}.summary > $dir/output/tts${counter}.ref.sorted.summary

        $pato $opts $flags -ds $tts_file -o $dir/output/tts${counter}.seq 1> /dev/null
        sort $dir/output/tts${counter}.seq.out > $dir/output/tts${counter}.seq.sorted.out
        sort $dir/output/tts${counter}.seq.summary > $dir/output/tts${counter}.seq.sorted.summary

        $pato_par $opts $flags -ds $tts_file -o $dir/output/tts${counter}.par 1> /dev/null
        sort $dir/output/tts${counter}.par.out > $dir/output/tts${counter}.par.sorted.out
        sort $dir/output/tts${counter}.par.summary > $dir/output/tts${counter}.par.sorted.summary

        diff $dir/output/tts${counter}.seq.sorted.out $dir/output/tts${counter}.ref.sorted.out > $dir/output/tts${counter}.seq.diff.out
        diff $dir/output/tts${counter}.par.sorted.out $dir/output/tts${counter}.ref.sorted.out > $dir/output/tts${counter}.par.diff.out
        diff $dir/output/tts${counter}.seq.sorted.summary $dir/output/tts${counter}.ref.sorted.summary > $dir/output/tts${counter}.seq.diff.summary
        diff $dir/output/tts${counter}.par.sorted.summary $dir/output/tts${counter}.ref.sorted.summary > $dir/output/tts${counter}.par.diff.summary

        result=0
        result=$(($result + $(wc -l < $dir/output/tts${counter}.seq.diff.out)))
        result=$(($result + $(wc -l < $dir/output/tts${counter}.par.diff.out)))
        result=$(($result + $(wc -l < $dir/output/tts${counter}.seq.diff.summary)))
        result=$(($result + $(wc -l < $dir/output/tts${counter}.par.diff.summary)))

        if [[ $result -eq 0 ]]; then
            printf "\r\e[1m[  \033[0;32mOK\033[0m  \e[1m]\e[0m\n"
        else
            printf "\r\e[1m[ \033[0;31mFAIL\033[0m \e[1m]\e[0m\n"
            failed=$(($failed + 1))
            fails=$((fails + 1))
        fi

        counter=$(($counter + 1))
    done
done
total=$(($total + $counter))
printf "         \e[1mRun:\e[0m ${counter} \e[1mOK\e[0m: $(($counter - $fails)) \e[1mFAIL:\e[0m ${fails}\n\n"

# TPX tests
declare -a tpx_opts=("" "-er 0" "-er 1" "-er 2" "-dd 1 -dc 2" "-e 10 -dd 2 -dc 3") 
declare -a tpx_flags=("" "-po true " "-of 1 " "-po true -of 1 " "-fr false " "-ssd false ")

fails=0
counter=0
printf "         \e[1mTPX search\e[0m\n"
for opts in "${tpx_opts[@]}"; do
    for flags in "${tpx_flags[@]}"; do
        printf "\e[1m[  --  ]\e[0m Test #${counter} ${flags}${opts}"

        sort $dir/ref/tpx${counter}.out > $dir/output/tpx${counter}.ref.sorted.out
        sort $dir/ref/tpx${counter}.summary > $dir/output/tpx${counter}.ref.sorted.summary

        $pato $opts $flags -ss $tfo_file -ds $tts_file -o $dir/output/tpx${counter}.seq 1> /dev/null
        sort $dir/output/tpx${counter}.seq.out > $dir/output/tpx${counter}.seq.sorted.out
        sort $dir/output/tpx${counter}.seq.summary > $dir/output/tpx${counter}.seq.sorted.summary

        $pato_par $opts $flags -ss $tfo_file -ds $tts_file -o $dir/output/tpx${counter}.par 1> /dev/null
        sort $dir/output/tpx${counter}.par.out > $dir/output/tpx${counter}.par.sorted.out
        sort $dir/output/tpx${counter}.par.summary > $dir/output/tpx${counter}.par.sorted.summary

        diff $dir/output/tpx${counter}.seq.sorted.out $dir/output/tpx${counter}.ref.sorted.out > $dir/output/tpx${counter}.seq.diff.out
        diff $dir/output/tpx${counter}.par.sorted.out $dir/output/tpx${counter}.ref.sorted.out > $dir/output/tpx${counter}.par.diff.out
        diff $dir/output/tpx${counter}.seq.sorted.summary $dir/output/tpx${counter}.ref.sorted.summary > $dir/output/tpx${counter}.seq.diff.summary
        diff $dir/output/tpx${counter}.par.sorted.summary $dir/output/tpx${counter}.ref.sorted.summary > $dir/output/tpx${counter}.par.diff.summary

        result=0
        result=$(($result + $(wc -l < $dir/output/tpx${counter}.seq.diff.out)))
        result=$(($result + $(wc -l < $dir/output/tpx${counter}.par.diff.out)))
        result=$(($result + $(wc -l < $dir/output/tpx${counter}.seq.diff.summary)))
        result=$(($result + $(wc -l < $dir/output/tpx${counter}.par.diff.summary)))

        if [[ $result -eq 0 ]]; then
            printf "\r\e[1m[  \033[0;32mOK\033[0m  \e[1m]\e[0m\n"
        else
            printf "\r\e[1m[ \033[0;31mFAIL\033[0m \e[1m]\e[0m\n"
            failed=$(($failed + 1))
            fails=$(($fails + 1))
        fi

        counter=$(($counter + 1))
    done
done
total=$(($total + $counter))
printf "         \e[1mRun:\e[0m ${counter} \e[1mOK\e[0m: $(($counter - $fails)) \e[1mFAIL:\e[0m ${fails}\n\n"

printf "\e[1mRun:\e[0m ${total} \e[1mOK\e[0m: $(($total - $failed)) \e[1mFAIL:\e[0m ${failed}\n"

# Clean up
rm -rf $dir/output

# Check if failed
if [[ $failed -gt 0 ]]
then
    exit 1
fi
exit 0
