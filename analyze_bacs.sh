#!/bin/bash

ASM=$1
BACS=bacs.fasta

# assert alignment
if [[ ! -f "$ASM.paf" ]]; then
  minimap2 --secondary=no -t 16 -ax asm20 $ASM $BACS > $ASM.sam
  samToErrorRate $ASM.sam $ASM | awk '{if ($9 == 0) print $0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t1\t"$12-$11"\t"$12-$10"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17}' > $ASM.paf
fi

# prep for calculating median, mean
cat <<EOF >Rstats
#!/usr/bin/R
dat = read.table("$ASM.tmp")
cat("Median:\t\t",  median(dat[,1]), "\n")
cat("MedianQV:\t", -10*log10(1-(median(dat[,1])/100)), "\n")
cat("Mean:\t\t", mean(dat[,1]), "\n")
cat("MeanQV:\t\t", -10*log10(1-(mean(dat[,1])/100)), "\n")
EOF

# get total BACs and BAC nucleotide
total=`grep -c ">" $BACS`
bp=`cat $BACS | grep --invert-match "^>" | wc -c`

# summarize BACs
echo "******************* BAC SUMMARY ******************" | tee $ASM.summary
echo " TOTAL    : $total" | tee -a $ASM.summary
echo " BP       : $bp" | tee -a $ASM.summary

# summarize assembly
echo "************** Processing $ASM.paf ****************"  | tee -a $ASM.summary

# get attempted BACs
cat $ASM.paf | awk '{if ($3 < -5000 && $4 > 90) print $1" "$2}' | sort | uniq | awk '{print $1}' | sort | uniq -c | awk '{if ($1 == 1) print $0}' >$ASM.attempted
ATMPT=`cat $ASM.attempted | wc -l`
cat $ASM.paf | awk '{if ($3 < -5000 && $4 > 90) print $1" "$2}' | sort | uniq | awk '{print $1}' | sort | uniq -c | awk '{if ($1 > 1) print $0}' >$ASM.multialign
MULTI=`cat $ASM.multialign | wc -l`

# get resolved BACs
cat $ASM.paf | awk '{if (($7-$6)/$8 > 0.995) print $1}'      | sort | uniq | awk '{print $1}' > $ASM.resolved
DONE=`cat $ASM.resolved | wc -l`
DONEBP=`cat $ASM.paf | awk '{if (($7-$6)/$8 > 0.995) print $1" "$8}' | sort | uniq | awk '{SUM+=$NF; } END { print SUM}'`

# get accuracies of resolved BACs
cat $ASM.paf | awk '{if (($7-$6)/$8 > 0.995) { print $0}}'  | awk '{print $4}' > $ASM.tmp
Rscript Rstats >$ASM.stats
rm Rstats $ASM.tmp

# print results and save
echo "$ATMPT $total" | awk      '{print "BACs attempted:     "$1" ("$1/$2")"} '  | tee -a $ASM.summary
echo "$MULTI $total $ATMPT" | awk      '{print "BACs multi-aligned: "$1" ("$1/$2") all aligned: "$1+$3" ("($1+$3)/$2")"} '  | tee -a $ASM.summary
echo "$DONE $total $DONEBP $bp $ATMPT" | awk '{print "BACs closed:        "$1" ("$1/$2" total, "$1/$5" attempted) BASES "$3" ("$3/$4" total)"}'  | tee -a $ASM.summary
cat $ASM.stats | tee -a $ASM.summary
echo "**********************************************"  | tee -a $ASM.summary

# save to table
echo -e "$ASM\t$total\t$ATMPT\t$DONE\t$bp\t$DONEBP\t$(cat $ASM.stats | sed 's/^.*:\s*//' | tr '\n' '\t')" > $ASM.table 
