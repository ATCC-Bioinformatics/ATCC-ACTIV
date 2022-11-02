#########################################
#   REQUIREMENTS:
#       lofreq, minimap2, bamqc/qualimap, samtools, kraken2, Nanofilt, python3 or greater
#
#   AUTHOR:
#       David A Yarmosh Jr., Senior Bioinformatician, ATCC
#
########################################

usage() {
    echo "Usage: run appropriate pipeline for submission
    -d kraken2 database; full path
    -l nanopore reads .fastq.gz; full path
    -o Output directory
    -r Reference fasta; full path
    -t Number of threads to use
    and anything else for help." 1>&2
    exit 1
}

while getopts 'd:l:o:r:t' OPTION
do
  case "$OPTION" in
    t) THREADS=$OPTARG
      ;;
    r) reference=$OPTARG
      ;;
    l) LONG=$OPTARG
      ;;
    d) DB=$OPTARG
      ;;
    f) outdir=$OPTARG
      ;;
    ?) usage
      exit 1
      ;;
  esac
done
shift "$(($OPTIND -1))"


bam_path=$outdir/$(basename $LONG .fastq.gz).bam
consensus_path=$outdir/$(basename $bam_path .bam).consensus.fasta
log_path=$outdir/run.log

#taxonomic binning with kraken2
bin_reads(){
  outdir=$1
  fwd=$2
  threads=$3
  bin1=$4
  db=$5
  echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
  echo "$(date) Beginning read classicfication with kraken." | tee -a "$log_path" >&2
  echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
  # Classify reads with kraken
  #### update to directory on HPC
  kraken2 --db $db \
  --threads "$threads" \
  --output "$outdir"/pre_QC/kraken_out.fastq \
  --report "$outdir"/pre_QC/kraken_report.txt \
  --gzip-compressed \
  $fwd
  if [ -f "$outdir"/pre_QC/kraken_out.fastq ]
  then
    echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
    echo "$(date) Read classification with kraken was successful." | tee -a "$log_path" >&2
    echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
  else
    echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
    echo "$(date) kraken output not found. Please check error log" | tee -a "$log_path" >&2
    echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
    exit
  fi

  # Bin on first bin taxid
  echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
  echo "$(date) Beginning read binning with extract_kraken_reads.py." | tee -a "$log_path" >&2
  echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2

  python /home/src/KrakenTools/extract_kraken_reads.py \
  -k "$outdir"/pre_QC/kraken_out.fastq \
  -r "$outdir"/pre_QC/kraken_report.txt \
  --fastq-output \
  -t $bin1 \
  --include-children \
  -s1 $fwd \
  -o  "$outdir"/$(basename $fwd .fastq.gz).binned.fastq
  if [ -f "$outdir"/$(basename $fwd .fastq.gz).binned.fastq ]
  then
    echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
    echo "$(date) Read binning with extract_kraken_reads.py was successful." | tee -a "$log_path" >&2
    echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
  else
    echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
    echo "$(date) extract_kraken_reads.py output not found. Please check error log" | tee -a "$log_path" >&2
    echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
    exit
  fi
  if [ ! -f "$outdir"/$(basename $fwd .fastq.gz).binned.fastq.gz ]
  then
  gzip "$outdir"/$(basename $fwd .fastq.gz).binned.fastq
  fi
}

map_reads() {
  local long_path=$1
  local reference_path=$2
  local output_dir_path=$3
  local threads=$4
  local bam_path=$5

  # variables for initial mapping and vcf used for bqsr
  local bam_pre_bqsr="$output_dir_path"/pre-bqsr.bam
  local bamqc_path="$output_dir_path"/bamqc
  local vcf_pre_bqsr="$output_dir_path"/pre-bqsr.vcf

  # variables for data generated by bqsr
  local cal_tbl_path="$output_dir_path"/bqsr_recalibration.tbl
  local cvtplot_path="$output_dir_path"/recalibration_quality.pdf # formerly, analyze_covariates

    echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
    echo "$(date) Starting dictionary creation for gatk." | tee -a "$log_path" >&2
    echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
    {
    # generate an initial bam file, sort it, assign a single read group, do qc
    # gatk requires explicit read groups for downstream steps
    #if the dict file exists, keep going. If it does not, run gatk CreateSequenceDictionary
    [ -f $(basename $reference_path .fasta).dict ] || gatk CreateSequenceDictionary -R "$reference_path"
    echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
    echo "$(date) Dictionary creation successful. Begin index creation and mapping with bwa mem." | tee -a "$log_path" >&2
    echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
    minimap2 -ax map-ont \
      -t $threads \
      "$reference_path" \
      "$long_path" |
    samtools sort -@ $threads - > $bam_path
    # these read group metadata are really just default placeholders
    }
      # ; } ||
    if [ -f $bam_path ]
    then
    echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
    echo "$(date) Mapping, realignment, sorting by leftmost coordinates, demultiplexing and artifact manipulation successful." | tee -a "$log_path" >&2
    echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
    else
    echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
    echo "$(date) Error during bwa index, bwa mem, lofreq viterbi, samtools sort, or picard.jar AddOrReplaceReadGroups. Quitting. Please check error log." |
    tee -a "$log_path" >&2
    echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
    exit
    fi

   echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
   echo "$(date) Starting initial variant calling + BQSR" | tee -a "$log_path" >&2
   echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
  { # bcftools without ploidy option for initial variant call set for BQSR
    qualimap bamqc \
      -nt "$threads" \
      -bam "$bam_path" \
      -outdir "$bamqc_path" \
      -outformat PDF:HTML
  }
      if [ -f "$bamqc_path"/report.pdf ]
      then
      echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
      echo "$(date) Initial variant calling, base quality score recalibration (BQSR), and statistical analysis of bam alignment successful." |
      tee -a "$log_path" >&2
      echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
      else
      echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
      echo "$(date) Error during initial variant calling (bcftools), BQSR (gatk), or bam quality control (qualimap bamqc). Quitting. Please check error log." |
      tee -a "$log_path" >&2
      echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
      exit
      fi
}

#call variants with freebayes
call_variants() {
  local reference_path=$1
  local bam_path=$2
  local outdir=$3

  mkdir -p $outdir/lofreq_filter
  name=$(basename $bam_path .bam)
  lofreq indelqual --dindel -f $reference_path $ban_path > $outdir/lofreq_filter/$name.indelqual.bam
  lofreq alnqual -b $outdir/lofreq_filter/$name.indelqual.bam $reference_path > $outdir/lofreq_filter/$name.indelqual.alnqual.bam
  lofreq call -f $reference_path --call-indels -C 50 $outdir/lofreq_filter/$name.indelqual.alnqual.bam > $outdir/lofreq_filter/$name.indelcall.vcf
  lofreq filter -i $outdir/lofreq_filter/$name.indelcall.vcf -v 50 -a 0.05 > $outdir/lofreq_filter/$name.final.vcf
}


#Now, the analysis begins
echo -e "$LONG\t$outdir"
longfq=$(basename $LONG)
gunzip -c $LONG | NanoFilt -q 15 -l 400 --headcrop 50 |gzip -c > "$outdir"/$longfq.filtered.fastq.gz #anticipates gzipped input files
bin_reads $outdir $outdir/$longfq.filtered.fastq.gz $THREADS 694009 $DB #SARS taxID, should contains SARS-CoV-2
medaka_haploid_variant -i $outdir/$longfq.filtered.binned.fastq.gz -r $reference -t $THREADS -o $outdir
call_variants $reference $bam_path $outdir
echo "Analysis complete at $outdir/lofreq_filter/$name.final.vcf"
