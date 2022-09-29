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
counter=0

# TFO tests
declare -a tfo_opts=("" "-l 14 -e 9" "-l 14 -e 10 -g 50 -G 50" "-l 14 -e 0 -m R" "-l 14 -e 0 -m Y" "-l 14 -e 0 -m M" "-e 10 -dl true -dd 1 -dc 2" "-e 10 -dl true -dd 2 -dc 2" "-e 10 -dl true -dd 1 -dc 2 -ssd off")
declare -a tfo_flags=("" "-po true " "-mf true ")

printf "Run mode: TFO search\n"
for opts in "${tfo_opts[@]}"; do
    for flags in "${tfo_flags[@]}"; do
        printf " [      ] Test #${counter} ${flags}${opts}"

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
            printf "\r [  \033[0;32mOK\033[0m  ]\n"
        else
            printf "\r [ \033[0;31mFAIL\033[0m ]\n"
            failed=$(($failed + 1))
        fi

        counter=$(($counter + 1))
    done
done

# Clean up
rm -rf $dir/output

# Check if failed
if [[ $FAILED -gt 0 ]]
then
    exit 1
fi
exit 0
