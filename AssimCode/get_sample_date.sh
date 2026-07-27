#!/bin/bash

read_month_day() {
    head -1 date_type.dat | awk '{print $2, $3}' | while read m d; do
        echo "${m} ${d}"
    done
}

IFS=' ' read current_month current_day < <(read_month_day)
declare -a month_days=(31 28 31 30 31 30 31 31 30 31 30 31)

generate_date_window() {
    local m=$(( 10#$1 )) d=$(( 10#$2 ))
    for offset in {-35..35..5}; do
        cm=$m cd=$d
        new_day=$((cd + offset))
        new_month=$cm
        
        while (( new_day > ${month_days[$new_month-1]} )); do
            new_day=$((new_day - ${month_days[$new_month-1]}))
            new_month=$((new_month % 12 + 1))
        done
        
        while (( new_day < 1 )); do
            new_month=$(( (new_month + 10) % 12 + 1 ))
            if (( new_month == 0 )); then
                new_month=12
            fi
            new_day=$((new_day + ${month_days[$new_month-1]}))
        done
        
        printf '%02d    %02d\n' "$new_month" "$new_day"
    done
}

> sample_date.dat
generate_date_window $current_month $current_day | sort -u -k1,1 -k2,2n >> sample_date.dat

