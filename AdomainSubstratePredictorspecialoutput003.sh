#!/bin/bash
#Use this script with caution, 'right' outcome is not guarranteed!!!
#This script should be used under a Unix-like system(eg, Ubuntu) with BLAST, HMMER (including Easel miniapps), bedtools, cd-hit, python3 and perl installed.
#You need to run 'chmod +x ./PKSclusterfinder.sh' to render the program executeble if this is its debut.
#intermediate--|-script
#              |-file
#Usage ./AdomainSubstratePredictor.sh PfamList(1) input_protein.faa(2) outputfolder(3) blastSeqDB(4) PfamKeyList(5)

#This is compiled by Chi Zhang @ Westlake University.

##################################
#00 quality control
##################################

#check if prameter number is right
if [ "$#" -ne 5 ]
then
  echo "Usage ./AdomainSubstratePredictor.sh PfamList(1) input_protein.faa(2) outputfolder(3) blastSeqDB(4) PfamKeyList(5)"
  exit 1
fi

:<<eof
#check if PfamList(1) meets Pfam.hmm format

if [ `grep -c "^HMMER" $1` -eq 0 -o `grep -c "^HMMER" $1` -ne `grep -c "^//" $1` ] 
then
  echo "Please check .hmm file!"
  exit 1
fi
eof

#check if input_protein.faa(2) meets fasta format

if [ `grep -c "^>" $2` -eq 0 ]
then
  echo "Please check .fa file!"
  exit 1
fi

#check the number of sequences in input_protein.faa(2)

echo "There are `grep -c "^>" $2` sequence(s) in submitted fasta file."


if [ `sed -e 's/^\(>[^[:space:]]*\).*/\1/' $2 | grep "^>" | sort -u | wc -l` -ne `sed -e 's/^\(>[^[:space:]]*\).*/\1/' $2 | grep -c "^>"` ]
then
  echo "Warning! There are duplicated ID(s) in fasta file and please remove duplication:"
  echo $(sed -e 's/^\(>[^[:space:]]*\).*/\1/' $2 | grep "^>" | sort | uniq -d)
  exit 1
fi

#check if no sequence has A domain at all

trap 'rm -f "$TMPFILE"' EXIT

TMPFILE=$(mktemp) || exit 1
hmmfetch $1 AMP-binding > $TMPFILE

if [ `hmmsearch --domE 1e-5 $TMPFILE $2 | grep -c "\[No hits detected that satisfy reporting thresholds\]"` -ne 0 ]
then
  echo "Please submit protein sequence(s) containing A domain(s)!"
  exit 1
fi

#check if outputfolder exist

if [ ! -e $3 ]
then
  echo "Outputfolder doesn't exist!."
  echo 'Do you want to create it (please make sure you have the privilege in parent directory if choose Y!)? (Y/N):'
  read choice
  if [ "$choice" = 'Y' -o "$choice" = 'y' ]
  then
    mkdir -p $3
    echo "$3 was created"
  else
    echo "Then, please specify another folder!"
    exit 1
  fi
fi

#create working subfolder and fetch matrix

WORKDIR=$3/NRPSPredictor$$ #remove last/ avoid directory error in python

mkdir -p $WORKDIR/00_QC/

sed -e 's/^\(>[^[:space:]]*\).*/\1/' $2 > $WORKDIR/00_QC/00_RAWFASTA.fa

sed -i -e 's/[[:space:]\/]/@/g' $WORKDIR/00_QC/00_RAWFASTA.fa

esl-sfetch --index $WORKDIR/00_QC/00_RAWFASTA.fa > /dev/null

##################################
#01 A domain (200AA) extraction and form a information table
##################################

mkdir $WORKDIR/01_DE/

hmmfetch $1 AdataSET0409200AA > $TMPFILE

hmmsearch --domT 40 --domtblout $WORKDIR/01_DE/200Ahmmsearch.dtbl --noali $TMPFILE $WORKDIR/00_QC/00_RAWFASTA.fa > /dev/null # change threshold 41 to 40

grep -v "^#" $WORKDIR/01_DE/200Ahmmsearch.dtbl | awk '{print $1"/"$20"-"$21, $20, $21, $1}' | sort -t'/' -k1,1 -k2n,2 | esl-sfetch -Cf $WORKDIR/00_QC/00_RAWFASTA.fa - | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' - | tail -n +2 > $WORKDIR/01_DE/ADomainhits.fa

grep ">" $WORKDIR/01_DE/ADomainhits.fa | awk -F'/' '{print $1}' - | sort -u > $WORKDIR/01_DE/ADomainhead.txt #just sequence ID not loc info

grep ">" $WORKDIR/01_DE/ADomainhits.fa > $WORKDIR/01_DE/ADomainHIThead.txt #seq ID with loc info

grep ">" $WORKDIR/00_QC/00_RAWFASTA.fa | sed -e 's/^\(>[^[:space:]]*\).*/\1/' | sort > $WORKDIR/01_DE/RAWhead.txt # all seq ID

comm -23 $WORKDIR/01_DE/RAWhead.txt $WORKDIR/01_DE/ADomainhead.txt > $WORKDIR/01_DE/NoADomainlist.txt

if [ `wc -l $WORKDIR/01_DE/NoADomainlist.txt | awk '{print $1}'` -ne 0 ] 
then
  echo "Warning! There is no A domain in the following sequences and they are not subject to upcoming analyses:"
  cat $WORKDIR/01_DE/NoADomainlist.txt
  echo
  echo
fi

##################################
#01A CD-HIT for time sake
################################## 

echo "CD-HIT MODE RUNNING ..."

cd-hit -i $WORKDIR/01_DE/ADomainhits.fa -o $WORKDIR/01_DE/ADomainhits0991.fa -c 0.991 -s 0.9 -T 0 -d 0 -M 0 > /dev/null # cd-hit for saving computational power 0.991 based on known threshold

echo `grep -c ">" $WORKDIR/01_DE/ADomainhits.fa`" A domains have been found and could be categorized into "`grep -c ">Cluster" $WORKDIR/01_DE/ADomainhits0991.fa.clstr`" clusters with 99.1% identity threshold."

echo "
inputfile1 = open('$WORKDIR/01_DE/ADomainhits0991.fa.clstr')
outputfile = open('$WORKDIR/01_DE/ADomainhits0991.fa.list', 'w')
line1 = inputfile1.readlines()
NUMLIST=[]
for i, element in enumerate(line1):
  if element.startswith('>'):
    NUMLIST.append(i)

xxx=0
LENGTH=len(NUMLIST)
while xxx < (LENGTH-1):
  if NUMLIST[xxx+1] - NUMLIST[xxx] > 1:
    REGION=[ [yyy.strip().split()[2][1:-3],yyy.strip().split()[-1]] for yyy in line1[NUMLIST[xxx]+1:NUMLIST[xxx+1]] ]
    xxx+=1
    for zzz in REGION:
      if zzz[1] == '*':
        AAA=zzz[0]
        REGION.remove(zzz)
        for UUU in REGION:
          outputfile.write(UUU[0]+'\t'+AAA+'\n')

" > $WORKDIR/01_DE/CDHITformatting.py

python3 $WORKDIR/01_DE/CDHITformatting.py

##################################
#02 A domain identity search (hmmalign+esl-alipid) IDP
################################## 

echo "IDP MODE RUNNING ..."

mkdir -p $WORKDIR/02_HMMALIGN/

# exit 1 # for test use

hmmalign --informat fasta --amino --outformat afa -o $WORKDIR/02_HMMALIGN/QUERYalignnoformat.fa $TMPFILE $WORKDIR/01_DE/ADomainhits.fa

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $WORKDIR/02_HMMALIGN/QUERYalignnoformat.fa | tail -n +2 > $WORKDIR/02_HMMALIGN/QUERYalign.fa

# cut -c2- $WORKDIR/01_DE/ADomainHIThead.txt > $WORKDIR/02_HMMALIGN/QUERYIDlist.txt

echo -e 'qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore' > $WORKDIR/02_HMMALIGN/TableHead

# blastp -query $WORKDIR/01_DE/ADomainhits.fa -db $4 -evalue 1e-5 -out $WORKDIR/02_HMMALIGN/IDPresult.txt -outfmt 6 -num_alignments 10 -num_threads 6

blastp -query $WORKDIR/01_DE/ADomainhits.fa -out $WORKDIR/02_HMMALIGN/BLASTresult.txt -task blastp-fast -db $4 -num_threads 12 -outfmt 6 -max_target_seqs 3 2> /dev/null

#grep -Fwf $WORKDIR/02_HMMALIGN/QUERYIDlist.txt $WORKDIR/02_HMMALIGN/IDPlist.txt | awk '{print $1,$2,$3}' > $WORKDIR/02_HMMALIGN/TargetIDPlist.txt

#cat $WORKDIR/02_HMMALIGN/TargetIDPlist.txt > $WORKDIR/02_HMMALIGN/IDPlist.txt # original IDPlist.txt is too large, so overwrite it with TargetIDPlist.txt after using


echo "IDP MODE DONE!"

##################################
#03 HMM whole sequence mode
##################################

echo "HMMwhole MODE RUNNING ..."

mkdir -p $WORKDIR/03_HMM/

echo -e 'qseqid\tHMMprofile\tBitscore\tlength\tE-value\tqstart\tqend' > $WORKDIR/03_HMM/TableHead

tail -n +1966 $1 > $WORKDIR/03_HMM/Onelayerforall.hmm

hmmsearch --domT 200 --domtblout $WORKDIR/03_HMM/NRPShmmsearch.dtbl --noali $WORKDIR/03_HMM/Onelayerforall.hmm $WORKDIR/01_DE/ADomainhits.fa > /dev/null

grep -v '^#' $WORKDIR/03_HMM/NRPShmmsearch.dtbl | awk '{print $1,$4,$8,$3,$7,$20,$21}' OFS='\t' | sort -k1,1 -k3nr,3 > $WORKDIR/03_HMM/ConciseHMMresult

echo "HMMwhole MODE DONE!"

##################################
#04 HMM key residues mode
##################################

echo "HMMkey MODE RUNNING ..."

mkdir -p $WORKDIR/04_KEYRES/

grep -v ">" $WORKDIR/02_HMMALIGN/QUERYalign.fa | cut -c7,19,21-23,25-39,46,48-66,82-84,100-104,123-128,130-133,136-137,139-140,149-153,155-156,159-163,177-180,183-189 > $WORKDIR/04_KEYRES/ONLYRES

sed -i -e "s%[-\.]%%g" $WORKDIR/04_KEYRES/ONLYRES #remove "-" in alignment

grep ">" $WORKDIR/02_HMMALIGN/QUERYalign.fa > $WORKDIR/04_KEYRES/ONLYHEAD

paste $WORKDIR/04_KEYRES/ONLYHEAD $WORKDIR/04_KEYRES/ONLYRES | tr "\t" "\n" > $WORKDIR/04_KEYRES/EXTRACTKEYRES.fa

hmmsearch --domtblout $WORKDIR/04_KEYRES/KEYREShmmsearch.dtbl --noali $5 $WORKDIR/04_KEYRES/EXTRACTKEYRES.fa > /dev/null

grep -v '^#' $WORKDIR/04_KEYRES/KEYREShmmsearch.dtbl | awk '{print $1,$4,$8,$3,$7,$20,$21}' OFS='\t' | sort -k1,1 -k3nr,3 > $WORKDIR/04_KEYRES/ConciseKEYRESresult

echo "HMMkey MODE DONE!"

##################################
#05 Outputformatting
##################################

echo "Last step: formatting..."

echo

mkdir -p $WORKDIR/05_Result/

# Here to compensate cd-hit results

cat $WORKDIR/02_HMMALIGN/TableHead $WORKDIR/02_HMMALIGN/BLASTresult.txt > $WORKDIR/05_Result/IDPresult.txt

cat $WORKDIR/03_HMM/TableHead $WORKDIR/03_HMM/ConciseHMMresult > $WORKDIR/05_Result/HMMwholeresult.txt

cat $WORKDIR/03_HMM/TableHead $WORKDIR/04_KEYRES/ConciseKEYRESresult > $WORKDIR/05_Result/HMMKEYresult.txt



mkdir $WORKDIR/05_Result/Formatting

touch $WORKDIR/05_Result/Formatting/UniqueIDlist.txt

grep ">" $WORKDIR/01_DE/ADomainhits.fa | sort -t'/' -k1,1 -k2n,2 -u | cut -c2- > $WORKDIR/05_Result/Formatting/UAHEADlist.txt

for UniqueID in $(cut -c2- $WORKDIR/01_DE/ADomainhead.txt)
do
	if [ `grep -w -c $UniqueID $WORKDIR/05_Result/Formatting/UniqueIDlist.txt` -eq 0 ]
	then
		echo $UniqueID >> $WORKDIR/05_Result/Formatting/UniqueIDlist.txt
	fi
done

echo -e 'ProteinID\t{HMMwhole}\t{HMMKEY}\t{IDP200}'

for ID in $(cat $WORKDIR/05_Result/Formatting/UniqueIDlist.txt)
do
  echo -e $ID'\t\c' >> $WORKDIR/05_Result/NRPSPredictorResult.txt
  grep -w $ID $WORKDIR/05_Result/Formatting/UAHEADlist.txt > $WORKDIR/05_Result/Formatting/${ID}AdomainHEAD.txt
  
  for subID in $(cat $WORKDIR/05_Result/Formatting/${ID}AdomainHEAD.txt)
  do
    # echo $subID >> $WORKDIR/05_Result/NRPSPredictorResult.txt

    # echo "[HMMwhole]" >> $WORKDIR/05_Result/NRPSPredictorResult.txt
    #substrate=`grep -w -m 1 $subID $WORKDIR/05_Result/HMMwholeresult.txt | cut -f2 | grep -o '{.*}'`
    if [ `grep -c -w $subID $WORKDIR/05_Result/HMMwholeresult.txt` -eq 0 ] 
    then
      echo -n '{unknown}' >> $WORKDIR/05_Result/NRPSPredictorResult.txt
    else
      echo -n `grep -w -m 1 $subID $WORKDIR/05_Result/HMMwholeresult.txt | awk '{print $2}' | sed 's/.*\({.*}\)/\1/'` >> $WORKDIR/05_Result/NRPSPredictorResult.txt
    fi
    
  done
  
  echo -e '\t\c' >> $WORKDIR/05_Result/NRPSPredictorResult.txt

  for subID in $(cat $WORKDIR/05_Result/Formatting/${ID}AdomainHEAD.txt)
  do

    # echo "[HMMkey]" >> $WORKDIR/05_Result/NRPSPredictorResult.txt
    if [ `grep -c -w $subID $WORKDIR/05_Result/HMMKEYresult.txt` -eq 0 ] 
    then
      echo -n '{unknown}' >> $WORKDIR/05_Result/NRPSPredictorResult.txt
    else
      echo -n `grep -w -m 1 $subID $WORKDIR/05_Result/HMMKEYresult.txt | awk '{print $2}' | sed 's/.*\({.*}\)/\1/'` >> $WORKDIR/05_Result/NRPSPredictorResult.txt
    fi
    
  done

  echo -e '\t\c' >> $WORKDIR/05_Result/NRPSPredictorResult.txt

  for subID in $(cat $WORKDIR/05_Result/Formatting/${ID}AdomainHEAD.txt)
  do
    # echo "[IDP]" >> $WORKDIR/05_Result/NRPSPredictorResult.txt
    echo -n `grep -w -m 1 $subID $WORKDIR/05_Result/IDPresult.txt | awk '{print $2}' | sed 's/.*\({.*}\)/\1/'` >> $WORKDIR/05_Result/NRPSPredictorResult.txt
  done
  
  echo -e '\n' >> $WORKDIR/05_Result/NRPSPredictorResult.txt
  echo >> $WORKDIR/05_Result/NRPSPredictorResult.txt

done

cat $WORKDIR/05_Result/NRPSPredictorResult.txt

echo "All done, cheer! Check result @ $WORKDIR/05_Result/NRPSPredictorResult.txt"

exit