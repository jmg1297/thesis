#!/bin/bash

# This should be run in the top-level projects directory

function fatalError {
echo $1
echo "Halting execution"
exit 1
}

if [ ! -f src/find_rmsk_gaps.py ];
then
fatalError "find_rmsk_gaps.py not found in $HOME/projects/src"
fi

RMSK=mm10_ref/repeatmasker/retrotransposons_only/rmskJoinedBaseline.RTs.txt
if [ ! -f $RMSK ];
then
fatalError $RMSK" not found"
fi

echo "Starting python script ... "
python src/find_rmsk_gaps.py $RMSK mm10_ref/repeatmasker/retrotransposon_indels/rt_indels.bed
echo "Done"

cd mm10_ref/repeatmasker/retrotransposon_indels

echo "Formatting rmsk names"

cat rt_indels.bed | \
sort -k1,1 -k2,2n | \
perl -lane 'if ($F[3] !~ /.*\/.*/){$F[3] =~ s/(.*)#(.*):(.*)/\1#\2\/\2:\3/} print join("\t", @F)' | \
sed -e 's/[#\/:]/ /g' -e 's/\?//g' | \
awk '{print $1,$2,$3,$5,$6,$4,$7,$8,$9}' | \
sed -r -e 's/ /:/4g' -e 's/:/ /4g' -e 's/\s+/\t/g' \
> rt_indels.formatted.bed

echo "Done"

echo "Extracting LINEs, SINEs, and LTRs individually"
grep -P "LINE" rt_indels.formatted.bed > LINE_indels.formatted.bed
grep -P "SINE" rt_indels.formatted.bed > SINE_indels.formatted.bed
grep -P "LTR" rt_indels.formatted.bed > LTR_indels.formatted.bed
echo "Done"
