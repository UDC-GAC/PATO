#!/usr/bin/env bash

# redirect all errors
exec 2> /dev/null

# set up
flavour=gnu
[ $# -gt 0 ] && flavour=$1

dir=$(dirname -- "$( readlink -f -- "$0"; )";)
rm -rf $dir/output
mkdir $dir/output

pato_seq=$dir/../target/$flavour/pato.serial
pato_par=$dir/../target/$flavour/pato.release

tfo_file=$dir/input/tfo.fa
tts_file=$dir/input/tts.fa

total=0
failed=0

# start the test suite
printf "\e[1mPATO: high PerformAnce TriplexatOr\e[0m -- Test suite\n\n"

# tfo tests
declare -a tfo_length_args=("" "-l 10 -L 14")
declare -a tfo_error_args=("" "-e 8 -E 10 -c 3 -g 50 -G 70")
declare -a tfo_match_args=("" "-b 3 -a on")
declare -a tfo_filter_args=("" "-fr off")
declare -a tfo_merge_args=("" "-mf on")
declare -a tfo_output_args=("" "-po on")

fails=0
counter=0
printf "         \e[1mTFO search\e[0m\n"
for len in "${tfo_length_args[@]}"; do
    for err in "${tfo_error_args[@]}"; do
        for mat in "${tfo_match_args[@]}"; do
            for fil in "${tfo_filter_args[@]}"; do
                for mer in "${tfo_merge_args[@]}"; do
                    for out in "${tfo_output_args[@]}"; do
                        args="${len} ${err} ${mat} ${fil} ${mer} ${out}"

                        printf "\e[1m[  --  ]\e[0m Test #${counter} ${args}"

                        sort $dir/ref/tfo${counter}.out > $dir/output/tfo${counter}.ref.sorted.out
                        sort $dir/ref/tfo${counter}.summary > $dir/output/tfo${counter}.ref.sorted.summary

                        $pato_seq $args -ss $tfo_file -o $dir/output/tfo${counter}.seq 1> /dev/null
                        sort $dir/output/tfo${counter}.seq.out > $dir/output/tfo${counter}.seq.sorted.out
                        sort $dir/output/tfo${counter}.seq.summary > $dir/output/tfo${counter}.seq.sorted.summary

                        $pato_par $args -ss $tfo_file -o $dir/output/tfo${counter}.par 1> /dev/null
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
                            rm -rf $dir/output/tfo${counter}*
                        else
                            printf "\r\e[1m[ \033[0;31mFAIL\033[0m \e[1m]\e[0m\n"
                            fails=$(($fails + 1))
                        fi

                        counter=$(($counter + 1))
                    done
                done
            done
        done
    done
done
printf "         \e[1mRun:\e[0m ${counter} \e[1mOK\e[0m: $(($counter - $fails)) \e[1mFAIL:\e[0m ${fails}\n\n"

total=$(($total + $counter))
failed=$(($failed + $fails))

# tts tests
declare -a tts_length_args=("" "-l 10 -L 14")
declare -a tts_error_args=("" "-e 8 -E 10 -c 3 -g 50 -G 70")
declare -a tts_match_args=("" "-b 3 -a on")
declare -a tts_filter_args=("" "-fr off")
declare -a tts_merge_args=("" "-mf on")
declare -a tts_output_args=("" "-po on")

fails=0
counter=0
printf "         \e[1mTTS search\e[0m\n"
for len in "${tts_length_args[@]}"; do
    for err in "${tts_error_args[@]}"; do
        for mat in "${tts_match_args[@]}"; do
            for fil in "${tts_filter_args[@]}"; do
                for mer in "${tts_merge_args[@]}"; do
                    for out in "${tts_output_args[@]}"; do
                        args="-cs 6 ${len} ${err} ${mat} ${fil} ${mer} ${out}"

                        printf "\e[1m[  --  ]\e[0m Test #${counter} ${args}"

                        sort $dir/ref/tts${counter}.out > $dir/output/tts${counter}.ref.sorted.out
                        sort $dir/ref/tts${counter}.summary > $dir/output/tts${counter}.ref.sorted.summary

                        $pato_seq $args -ds $tts_file -o $dir/output/tts${counter}.seq 1> /dev/null
                        sort $dir/output/tts${counter}.seq.out > $dir/output/tts${counter}.seq.sorted.out
                        sort $dir/output/tts${counter}.seq.summary > $dir/output/tts${counter}.seq.sorted.summary

                        $pato_par $args -ds $tts_file -o $dir/output/tts${counter}.par 1> /dev/null
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
                            rm -rf $dir/output/tts${counter}*
                        else
                            printf "\r\e[1m[ \033[0;31mFAIL\033[0m \e[1m]\e[0m\n"
                            fails=$(($fails + 1))
                        fi

                        counter=$(($counter + 1))
                    done
                done
            done
        done
    done
done
printf "         \e[1mRun:\e[0m ${counter} \e[1mOK\e[0m: $(($counter - $fails)) \e[1mFAIL:\e[0m ${fails}\n\n"

total=$(($total + $counter))
failed=$(($failed + $fails))

# tpx tests
declare -a tpx_length_args=("" "-l 10 -L 14")
declare -a tpx_error_args=("" "-e 8 -E 10 -c 3 -g 50 -G 70")
declare -a tpx_match_args=("" "-b 3 -a on")
declare -a tpx_filter_args=("" "-fr off")
declare -a tpx_merge_args=("" "-er 1" "-er 2")
declare -a tpx_output_args=("" "-of 1")

fails=0
counter=0
printf "         \e[1mTPX search\e[0m\n"
for len in "${tpx_length_args[@]}"; do
    for err in "${tpx_error_args[@]}"; do
        for mat in "${tpx_match_args[@]}"; do
            for fil in "${tpx_filter_args[@]}"; do
                for mer in "${tpx_merge_args[@]}"; do
                    for out in "${tpx_output_args[@]}"; do
                        args="-cs 4 ${len} ${err} ${mat} ${fil} ${mer} ${out}"

                        printf "\e[1m[  --  ]\e[0m Test #${counter} ${args}"

                        sort $dir/ref/tpx${counter}.out > $dir/output/tpx${counter}.ref.sorted.out
                        sort $dir/ref/tpx${counter}.summary > $dir/output/tpx${counter}.ref.sorted.summary

                        $pato_seq $args -ss $tfo_file -ds $tts_file -o $dir/output/tpx${counter}.seq 1> /dev/null
                        sort $dir/output/tpx${counter}.seq.out > $dir/output/tpx${counter}.seq.sorted.out
                        sort $dir/output/tpx${counter}.seq.summary > $dir/output/tpx${counter}.seq.sorted.summary

                        $pato_par $args -ss $tfo_file -ds $tts_file -o $dir/output/tpx${counter}.par 1> /dev/null
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
                            rm -rf $dir/output/tpx${counter}*
                        else
                            printf "\r\e[1m[ \033[0;31mFAIL\033[0m \e[1m]\e[0m\n"
                            fails=$(($fails + 1))
                        fi

                        counter=$(($counter + 1))
                    done
                done
            done
        done
    done
done
printf "         \e[1mRun:\e[0m ${counter} \e[1mOK\e[0m: $(($counter - $fails)) \e[1mFAIL:\e[0m ${fails}\n\n"

total=$(($total + $counter))
failed=$(($failed + $fails))

# end the test suite and report the results
printf "\e[1mRun:\e[0m ${total} \e[1mOK\e[0m: $(($total - $failed)) \e[1mFAIL:\e[0m ${failed}\n"

# did it fail
if [[ $failed -gt 0 ]]; then
    exit 1
fi

# clean up
rm -rf $dir/output

exit 0
