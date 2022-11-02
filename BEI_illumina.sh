#########################################
#	REQUIREMENTS:
#		freebayes, lofreq, bwa, gatk, bamqc/qualimap, bcftools, samtools, picardtools, kraken2, fastp, python3 or greater, seqkit, vt
#
#	AUTHOR:
#		David A Yarmosh Jr., Senior Bioinformatician, ATCC
#
########################################
usage() {
    echo "Usage: run appropriate pipeline for submission
    -1 Illumina forward reads .fastq.gz; full path
    -2 Illumina reverse reads .fastq.gz; full path
    -d kraken2 database; full path
    -o Output directory
    -r Reference fasta; full path
    -t Number of threads to use
    and anything else for help." 1>&2
    exit 1
}

while getopts '1:2:d:o:r:t' OPTION
do
  case "$OPTION" in
    t) THREADS=$OPTARG
      ;;
    r) reference=$OPTARG
      ;;
    1) FWD=$OPTARG
      ;;
    2) REV=$OPTARG
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

bam_path=$outdir/$(basename $FWD _1.fastq).bam
vcf_path=$outdir/$(basename $bam_path .bam).vcf
simplified_vcf_path=$vcf_path.simplified.vcf
consensus_path=$outdir/$(basename $bam_path .bam).consensus.fasta
log_path=$outdir/run.log



#pre-QC function
pre_QC(){
  outdir=$1
  threads=$2
  fwd=$3
  rev=$4
  mkdir -p $outdir/pre_QC
  echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
  echo "$(date) Beginning Pre-QC with fastp." | tee -a "$log_path" >&2
  echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
  #### Pre-QC ####
  # fastp --n_base_limit 10 --qualified_quality_phred 30 --unqualified_percent_limit 10 \
  # -i $fwd -I $rev -o "$outdir"/$(basename $fwd .fastq.gz).filtered.fastq.gz -O "$outdir"/$(basename $rev .fastq.gz).filtered.fastq.gz \
  # -h "$outdir"/pre_QC/fastp.html -j "$outdir"/pre_QC/fastp.json
  # fastp --n_base_limit 10 --qualified_quality_phred 25 --unqualified_percent_limit 20 \
  # -i $fwd -I $rev -o "$outdir"/$(basename $fwd .fastq.gz).filtered.fastq.gz -O "$outdir"/$(basename $rev .fastq.gz).filtered.fastq.gz \
  # -h "$outdir"/pre_QC/fastp.html -j "$outdir"/pre_QC/fastp.json
  fastp -i $fwd -I $rev -o "$outdir"/$(basename $fwd .fastq.gz).filtered.fastq.gz \
  --detect_adapter_for_pe \
  -O "$outdir"/$(basename $rev .fastq.gz).filtered.fastq.gz \
  -h "$outdir"/pre_QC/fastp.html -j "$outdir"/pre_QC/fastp.json
  if [ -f "$outdir"/$(basename $rev .fastq.gz).filtered.fastq.gz ]
  then
    echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
    echo "$(date) Pre-QC with fastp successful." | tee -a "$log_path" >&2
    echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
  else
    echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
    echo "$(date) fastp output not found. Please check error log" | tee -a "$log_path" >&2
    echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
    exit
  fi
}

#taxonomic binning with kraken2
#requires kraken2 installed, its bacterial_viral_db has a hardcoded location here that needs to be changed, and KrakenTools has a hardcoded location for extract_kraken_reads.py that also needs to be changed.
bin_reads(){
  outdir=$1
  fwd=$2
  rev=$3
  threads=$4
  bin1=$5
  db=$6
  echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
  echo "$(date) Beginning read classicfication with kraken." | tee -a "$log_path" >&2
  echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
  # Classify reads with kraken
  #### update to directory on HPC
  kraken2 --db /home/shared/databases/kraken/bacterial_viral_db \
  --threads "$threads" \
  --output "$outdir"/pre_QC/kraken_out.fastq \
  --report "$outdir"/pre_QC/kraken_report.txt \
  --gzip-compressed \
  --paired $fwd $rev
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
  #### update to directory on HPC
  python /home/src/KrakenTools/extract_kraken_reads.py \
  -k "$outdir"/pre_QC/kraken_out.fastq \
  -r "$outdir"/pre_QC/kraken_report.txt \
  --fastq-output \
  -t $bin1 \
  --include-children \
  -s1 $fwd \
  -s2 $rev \
  -o  "$outdir"/$(basename $fwd .fastq.gz).binned.fastq \
  -o2 "$outdir"/$(basename $rev .fastq.gz).binned.fastq
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
  gzip "$outdir"/$(basename $rev .fastq.gz).binned.fastq
#no guarantee kraken will extract the same reads from forward and reverse, so seqkit common is used to ensure identical reads are present in both files.
  seqkit common $outdir/$fwd2.filtered.binned.fastq.gz $outdir/$rev2.filtered.binned.fastq.gz | gzip > $outdir/$fwd2.filtered.binned.common.fastq.gz
  seqkit common $outdir/$rev2.filtered.binned.fastq.gz $outdir/$fwd2.filtered.binned.fastq.gz | gzip > $outdir/$rev2.filtered.binned.common.fastq.gz

  fi
}

#map reads to reference, realign with lofreq, and recalibrate quality scores
map_reads() {
  local r1_path=$1
  local r2_path=$2
  local reference_path=$3
  local output_dir_path=$4
  local threads=$5
  local bam_path=$6

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
    bwa index "$reference_path" &&
    bwa mem \
      -t $threads \
      "$reference_path" \
      "$r1_path" \
      "$r2_path" |
    lofreq viterbi -f "$reference_path" - |
    samtools sort -@ $threads - |
    # these read group metadata are really just default placeholders
    java -jar /home/src/picard.jar AddOrReplaceReadGroups \
      -I /dev/stdin \
      -O "$bam_pre_bqsr" \
      -RGLB 1 \
      -RGPL illumina \
      -RGPU unit1 \
      -RGSM sample
    }
      # ; } ||
    if [ -f $bam_pre_bqsr ]
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
    bcftools mpileup \
      -f "$reference_path" \
      -d 8000 \
      --threads "$threads" \
      -Ou \
      "$bam_pre_bqsr" |
    bcftools call --threads "$threads" -M -mv -Ou |
    bcftools filter --threads "$threads" -Oz -i '%QUAL>=30' -o "$vcf_pre_bqsr".gz &&
    /home/src/gatk-4.2.2.0/gatk IndexFeatureFile -I "$vcf_pre_bqsr".gz --verbosity DEBUG &&
    # GATK run for BQSR
    /home/src/gatk-4.2.2.0/gatk BaseRecalibrator \
      -I "$bam_pre_bqsr" \
      -R "$reference_path" \
      --known-sites "$vcf_pre_bqsr".gz \
      --verbosity DEBUG \
      -O "$cal_tbl_path" &&
    /home/src/gatk-4.2.2.0/gatk ApplyBQSR \
      -R "$reference_path" \
      -I "$bam_pre_bqsr" \
      --bqsr-recal-file "$cal_tbl_path" \
      --verbosity DEBUG \
      -O "$bam_path" &&
    /home/src/gatk-4.2.2.0/gatk AnalyzeCovariates \
      -bqsr "$cal_tbl_path" \
      --verbosity DEBUG \
      -plots "$cvtplot_path" &&
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

#variant calling with lofreq
call_variants() {
  local reference_path=$1
  local bam_path=$2
  local vcf_path=$3
  local output_dir=$4
  mkdir -p $output_dir/lofreq_filter
  name=$(basename $bam_path .bam)
  lofreq indelqual --dindel -f $reference_path $bam_path > $output_dir/lofreq_filter/$name.indelqual.bam
  lofreq alnqual -b $output_dir/lofreq_filter/$name.indelqual.bam $reference_path > $output_dir/lofreq_filter/$name.indelqual.alnqual.bam
  lofreq call -f $reference_path --call-indels -C 50 $output_dir/lofreq_filter/$name.indelqual.alnqual.bam > $output_dir/lofreq_filter/$name.indelcall.vcf
  lofreq filter -i $output_dir/lofreq_filter/$name.indelcall.vcf -v 50 -a 0.15 > $vcf_path

  gzip -c "$vcf_path" > "$vcf_path".gz
}

#Analysis begins here:
echo -e "$FWD\t$REV\t$outdir"

pre_QC $outdir $THREADS $FWD $REV
fwd2=$(basename $FWD .fastq.gz)
rev2=$(basename $REV .fastq.gz)
bin_reads $outdir $outdir/$fwd2.filtered.fastq.gz $outdir/$rev2.filtered.fastq.gz $THREADS 694009 $DB #SARS taxID, should contains SARS-CoV-2
map_reads $outdir/$fwd2.filtered.binned.common.fastq.gz $outdir/$rev2.filtered.binned.common.fastq.gz $reference $outdir $THREADS $bam_path
call_variants $reference $bam_path $vcf_path $outdir
