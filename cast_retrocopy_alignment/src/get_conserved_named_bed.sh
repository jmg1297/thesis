#!/bin/bash

INPUT=$1
OUTPUT=$2

cat $INPUT | \
awk '{print $1,$2,$3,$6,0,$4}' | \
perl -lane 'BEGIN{my %retro; my $count; my $name} $name=$F[3]; if (defined $retro{$name}){$retro{$name} += 1; $count=$retro{$name}} else {$count=1; $retro{$name} = 1}; $F[3]="retro-$name-$count"; print "@F"' | \
sed -r 's/\s+/\t/g' \
> $OUTPUT
