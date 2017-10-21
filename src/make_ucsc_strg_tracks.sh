#!/bin/bash

function gtfToBigBed {
	gtfToGenePred -genePredExt $1/out.gtf $1/out.genePred.chr;
	perl -lane 'if ($F[1] !~ /^[(GL)|(JH)|(MT)]/) {$F[1] = "chr".$F[1]; print join "\t", @F}' $1/out.genePred.chr > $1/out.genePred
	rm $1/out.genePred.chr
	genePredToBigGenePred $1/out.genePred stdout | sort -k1,1 -k2,2n > $1/out.bigGenePred
	bedToBigBed -type=bed12+8 -tab -as=$2 $1/out.bigGenePred $3 $4/$5.bigBed
}

function appendToTrackDb {
	SAMPLE=$2
	echo "track "$SAMPLE >> $1
	echo "bigDataUrl http://jeremy-bio.sysbiol.cam.ac.uk/JG600_UCSC/TRANSCRIPTOMES/mm10/"$SAMPLE".bigBed" >> $1
	echo "shortLabel "$SAMPLE >> $1
	echo "longLabel "$SAMPLE" StringTie transcriptome" >> $1
	echo "type bigGenePred" >> $1
	echo "" >> $1
}

#gtfToBigGenePred /home/jg600/projects/rna-seq_transcriptome_reconstruction/data/stringtie/10966_8#1_C57BL-6J_Female_B /home/jg600/projects/rna-seq_transcriptome_reconstruction/data/stringtie/10966_8#1_C57BL-6J_Female_B/bigGenePred.as /home/jg600/projects/mm10_ref/genome/mm10.chrom.sizes

for dir in $(find $1 -type d | sed '1d');
do
	SAMPLE=$(echo $dir | awk 'BEGIN{FS="/"} {print $NF}' | sed -r 's/#/-/g')
	echo $SAMPLE
	gtfToBigBed $dir $HOME/projects/mm10_ref/bigGenePred.as $HOME/projects/mm10_ref/genome/mm10.chrom.sizes $HOME/UCSC_HUBS/TRANSCRIPTOMES/mm10 $SAMPLE
	appendToTrackDb $HOME/UCSC_HUBS/TRANSCRIPTOMES/mm10/trackDb.txt $SAMPLE 
done
