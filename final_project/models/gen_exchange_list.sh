#!/bin/bash

FN="-"

if [ -n "$1" ]
then
	FN="$1"
fi

cat "$FN" | grep "EX_" | sed 's/^[[:blank:]]*<reaction //g; s/>$//g; s/fbc://g' | while read line; do eval $line; echo -e "$id\t $name"; done
