#!/usr/bin/env bash

MMGENOME=/nfs1/MICRO/Dreher_Lab/programs/mmgenome

show_usage () { 
	echo -e "\nThis script runs the data generation workflow for mmgenome.\n\nThe script can be used in two modes:\n\t(1) full: -m f\n\t(2) coverage only: -m c\n\nPlease feed files to the script in the order shown below.\n\nUsage:\nmm.datagen.sh -m [f|c] -t [number of threads (integer)] -n [name of reads] <Assembly file (fasta)> <Forward read file (fastq)> <Reverse read file (fastq)>\n\nExample:\nmm.datagen.sh -m f -t 12 -n wa102c wa102c.fa R1.fastq R2.fastq\n\nWARNING: DO NOT USE hyphens \"-\" or pluses \"+\" in the name of your assembly file.  Phylopythia will not work if you do . . . and it will take you hours to find out.\n\n"
}

if [ "$#" == 0 ]; then
	show_usage
	exit 1
fi

clear

echo "---Metagenomics workflow script v.2.1.1 for Dreher Lab---"

while getopts ":m:t:n:" OPTS; do
	case $OPTS in
		m)
			if [[ $OPTARG =~ [Ff] ]]; then
				FULL=1
				echo -e "\nRunning in \"full output\" mode; generating all output files."  >&2
			elif [[ $OPTARG =~ [Cc] ]]; then
				FULL=0
				echo -e "\nRunning in \"coverage only\" mode; generating only coverage output file."  >&2
			elif [[ $OPTARG =~ [FfCc] ]]; then
				echo -e "Invalid option: -$OPTARG\n"
			else
				break
			fi
			;;&
		t)
			if [[ $OPTARG =~ [[:digit:]] ]]; then			
				THREADS=$OPTARG
			else [[ "$OPTARG" =~ [[:digit:]] ]]
				echo -e "Invalid option: -$OPTARG\n"
			fi
			;;&
                n)
                  	if [[ "$OPTARG" =~ [[:alnum:]] ]]; then
                                        READS=$OPTARG
                        else
                            	echo -e "Invalid argument: -n requires a name for the reads" >&2
                        fi
                        ;;&
		\?)
			echo -e "Invalid option: -$OPTARG\n" >&2			
			show_usage
			exit 2
			;;
		:)
			echo -e "Option -$OPTARG needs an argument.\n" >&2
			show_usage
			exit 3
			;;
	esac
done

STEM=${7%.*}
RENAME=$STEM.rn.fa
R1=$8
R2=$9
OUTPUT="$READS.$STEM"

mkdir output 2> /dev/null
mkdir working 2> /dev/null

echo -e "\nRenaming each of the contigs in the assembly with a unique number"
awk '/>/{print ">"i++; next} {print}' $7 > output/$RENAME
cp  output/$RENAME working/$RENAME
echo -e "\nWriting a PhylopythiaS+ configuration file and execution script ("ppsp.sh") in a sub-directory for YOU to run later ;)"
mkdir ppsp_$STEM 2> /dev/null
echo -e "[PhyloPythiaS_Plus]\npipelineDir=$PWD/ppsp_$STEM\ninputFastaFile=$PWD/output/$RENAME" > ppsp_$STEM/ppsp_$STEM.cfg
echo -e "inputFastaScaffoldsFile=\nscaffoldsToContigsMapFile=\nreferencePlacementFileOut=\ndatabaseFile=/nfs1/MICRO/Dreher_Lab/db/ppsp/reference_NCBI20140513/taxonomy\nrefSeq=/nfs1/MICRO/Dreher_Lab/db/ppsp/reference_NCBI20140513/centroids_noplasmids\nexcludeRefSeqRank=\ns16Database=/nfs1/MICRO/Dreher_Lab/db/ppsp/reference_NCBI20140513/silva115\nmgDatabase=/nfs1/MICRO/Dreher_Lab/db/ppsp/reference_NCBI20140513/mg4\nexcludeRefMgRank=\nppsInstallDir=/nfs1/MICRO/Dreher_Lab/programs/pps/tools_1427936977/PhyloPythiaS/1_3\nmothurInstallDir=/nfs1/MICRO/Dreher_Lab/programs/pps/tools_1427936977/mothur/mothur1333\nhmmerBinDir=/nfs1/MICRO/Dreher_Lab/programs/pps/tools_1427936977/hmmer-3.0/binaries\nrnaHmmInstallDir=/nfs1/MICRO/Dreher_Lab/programs/pps/tools_1427936977/rna_hmm3\nrankIdCut=6\nmaxLeafClades=50\nminPercentInLeaf=1.0\nminSeqLen=1000\nplaceContigsFromTheSameScaffold=True\nagThreshold=0.3\nassignedPartThreshold=0.5\nrankIdAll=0\nweightStayAll=60.0\nrankIdCutMinBp=10000\noverrideWithPPSPlacements=True\nminBpToModel=100000\nminSSDfileSize=5000\nmaxSSDfileSize=400000\nminGenomesWgs=1\nminBpPerSpecies=300000\ncandidatePlTopPercentThreshold=0.1\noutputFileContigSubPattern=^(.*)\nmothurClassifyParamOther=method=bayesian, cutoff=80, iters=300\nconfigPPS=\nrecallMinFracClade=0.001\nprecisionMinFracPred=0.001\ncorrectlabelthreshold=0.9" >> ppsp_$STEM/ppsp_$STEM.cfg
echo -e "#!/usr/bin/env bash\nppsp -c ppsp_$STEM.cfg -n -g -o s16 mg -t -p c -r -s\ncp output/$RENAME.PP.pOUT ../output/pps.txt" > ppsp_$STEM/ppsp.sh
chmod +x ppsp_$STEM/ppsp.sh

echo -e "\nMapping reads to contigs"
BAM=$OUTPUT.bam
SAM=$OUTPUT.sam

bwa index working/$RENAME 2> /dev/null
bwa mem -t $THREADS working/$RENAME $R1 $R2 2> /dev/null | samtools view -bS - | samtools sort - working/$OUTPUT
samtools index working/$BAM
samtools view working/$BAM > working/$SAM

echo -e "\nDetermining average read coverage of contigs"
samtools depth working/$BAM | $MMGENOME/scripts/calc.coverage.in.bam.depth.pl -i - -o output/$OUTPUT.cov

if [ "$FULL" == 0 ]; then
	exit 6
fi

echo -e "\nWriting out paired-end network connections between contigs if they have three or more connections."
$MMGENOME/scripts/network.pl -i working/$SAM -o output/network.txt -f 3 2> /dev/null

echo -e "\nPredicting proteins in the assembly contigs"
prodigal -a working/$STEM.temp.orfs.faa -i output/$RENAME -m -o working/$STEM.temp.txt -p meta -q
cut -f1 -d " " working/$STEM.temp.orfs.faa > working/$STEM.orfs.faa

echo -e "\nFinding essential genes using an HMM search"
hmmsearch --cpu $THREADS --tblout working/$STEM.hmm.orfs.txt --cut_tc --notextw $MMGENOME/scripts/essential.hmm working/$STEM.orfs.faa > working/$STEM.hmm.temp.txt
echo "scaffold orf hmm.id" > output/essential.txt
tail -n+4  working/$STEM.hmm.orfs.txt | sed 's/ * / /g' | cut -f1,4 -d " " | sed 's/_/ /' >> output/essential.txt
grep -v "#" working/$STEM.hmm.orfs.txt | cut -f1 -d " " > working/$STEM.list.of.positive.orfs.txt
perl $MMGENOME/scripts/extract.using.header.list.pl -l working/$STEM.list.of.positive.orfs.txt -s working/$STEM.orfs.faa -o working/$STEM.orfs.hmm.faa

echo -e "\nIdentifying essential genes using the refseq protein database"
#blastp -query working/$STEM.orfs.hmm.faa -db refseq_protein -evalue 1e-5 -num_threads $THREADS -max_target_seqs 5 -outfmt 5 -out working/$STEM.orfs.hmm.blast.xml
rapsearch -q working/$STEM.orfs.hmm.faa -d /nfs1/MICRO/Dreher_Lab/db/rap2_refseq -o working/$STEM.orfs.hmm.rapsearch -z $THREADS -v 5 -b 5 -t a -x t

echo -e "\nExtracting the consensus taxonomic classification of the essential genes"
echo -e "load taxGIFile='/nfs1/MICRO/Dreher_Lab/programs/megan5/gi_taxid-March2015X.bin'\nimport blastfile='working/$STEM.orfs.hmm.rapsearch.xml' meganfile='working/temp.rma' blastFormat=BlastXML;\nrecompute toppercent=5;\nrecompute minsupport=1;\nupdate;\ncollapse rank=Species;\nupdate;\nselect nodes=all;\nexport what=CSV format=readname_taxonpath separator=tab file='working/$STEM.orfs.hmm.rapsearch.tax.txt';\nupdate;\nclose" > working/megan_commands.txt
MEGAN5 < working/megan_commands.txt
perl $MMGENOME/scripts/hmm.majority.vote.pl -i working/$STEM.orfs.hmm.rapsearch.tax.txt -o output/tax.txt
