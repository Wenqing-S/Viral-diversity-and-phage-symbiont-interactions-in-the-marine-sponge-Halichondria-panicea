######################
### short reads, contigs and MAGs
metawrap read_qc -1 00_RAW_DATA/read_1.fastq -2 00_RAW_DATA/read_2.fastq -t 32 -o 01_READ_QC -x HP_genome
 mv 01_READ_QC/final_pure_reads_1.fastq 02_CLEAN_READS/read_1.fastq 
 mv 01_READ_QC/final_pure_reads_2.fastq 02_CLEAN_READS/read_2.fastq
metawrap assembly \
   -1 02_CLEAN_READS/read_1.fastq \
   -2 02_CLEAN_READS/read_2.fastq \
   -m 800 -t 32  \
   -o 03_ASSEMBLY
metawrap binning \
  -o 05_INITIAL_BINNING \
  -t 32 \
	  -a /03_ASSEMBLY/megahit_coassembly/final.contigs.fa \
  --metabat2 \
  --maxbin2 \
  --concoct \
  02_CLEAN_READS/*.fastq
metawrap bin_refinement \
  -o 06_BIN_REFINEMENT \
  -t 32 \
  -A 05_INITIAL_BINNING/metabat2_bins/ \
  -B 05_INITIAL_BINNING/maxbin2_bins/ \
  -C 05_INITIAL_BINNING/concoct_bins/ \
  -c 50 \
  -x 10
metawrap reassemble_bins \
  -o 09_BIN_REASSEMBLY \
  -1 $r1_files \
  -2 $r2_files \
  -t 20 \
  -m 800 \
  -c 50 \
  -x 10 \
  -b 06_BIN_REFINEMENT/metawrap_50_10_bins
######################
### long reads, contigs and MAGs
  hifiasm_meta -t 32 -o asm_output ERR11725114.fastq.gz
  flye -o PacBio/01_assembly/metaflye --meta --pacbio-hifi ERR11725114.fastq.gz --keep-haplotypes --hifi-error 0.01 -t 32
  python /MetaCC/MetaCC.py norm -e DpnII -e HinfI -v \
    /PacBio/01_assembly/hifiasm-meta/hifiasm.fasta \
    /PacBio/03_binning/HiCBin/hp/hic_mapping_output/MAP_89_SORTED.bam \
    /PacBio/03_binning/MetaCC/89_OUTPUT
  python /MetaCC/MetaCC.py bin --cover -v /PacBio/01_assembly/hifiasm-meta/hifiasm.fasta /PacBio/03_binning/MetaCC/89_OUTPUT
#also use HiFi-MAG-Pipeline get bins from metabat2&semibin2
metawrap bin_refinement \
  -o 06_BIN_REFINEMENT \
  -t 32 \
  -A metabat2_bins/ \
  -B metacc_bins/ \
  -C semi_bins/ \
  -c 50 \
  -x 10
######################
### quality and taxonomy of MAGs 
dRep dereplicate MAG_drep -g /MAG/*.fa -p 32 -comp 50 -con 10
checkm2 predict --threads 30 --input BIN --output-directory checkm2
gtdbtk classify_wf \
  --genome_dir metawrap_50_10_bins \
  --extension fa \
  --out_dir gtdbtk_output \
  --cpus 12
######################
### microbial diversity
phyloFlash.pl -dbhome /database/138.1 -lib 373 -CPUs 16 \
 -read1 /02_CLEAN_READS/read_1.fastq \
 -read2 /02_CLEAN_READS/read_2.fastq \
 -almosteverything
iqtree2 \
  -s bac120_user_msa.fasta \
  -m MFP \
  -B 1000 \
  -T 16
  
######################
### viral analysis
### identification
ViWrap run --input_metagenome /03_ASSEMBLY/final_assembly.fasta \
           --input_reads /02_CLEAN_READS/read_1.fastq,02_CLEAN_READS/read_2.fastq \
           --out_dir /viwrap \
           --db_dir /database/viwrap/ViWrap_db \
           --identify_method vb-vs-dvf \
           --threads 32 \
           --input_length_limit 5000       \
           --conda_env_dir /database/viwrap/env

makeblastdb -in /All_virus/all_virus.fasta -dbtype nucl -out votu
blastn -query All_virus/all_virus.fasta -db votu_all -outfmt '6 std qlen slen' -max_target_seqs 10000 -out my_blast.tsv -num_threads 32
/checkv_script/anicalc.py -i my_blast.tsv -o my_ani.tsv
/checkv_script/aniclust.py --fna all_virus.fasta --ani my_ani.tsv --out my_clusters.tsv --min_ani 95 --min_tcov 85 --min_qcov 0
seqtk subseq All_virus/all_virus.fasta votu_name.txt > vOTU.fasta
checkv end_to_end vOTU.fasta vOTU_checkv -d /CheckV_db -t 16
genomad end-to-end --cleanup --splits 32 vOTU.fasta vOTU_genomad_output /database/genomad/genomad_db
vcontact3 run --nucleotide /vOTU.fasta --output vcontact3 --db-version 228 --db-path /database/vcontact3
### annotation
virsorter run --prep-for-dramv -w vs2.out -i vOTU.fasta -j 32 all
DRAM-v.py annotate -i for-dramv/final-viral-combined-for-dramv.fa \
       -v /for-dramv/viral-affi-contigs-for-dramv.tab -o dramv_result --threads 16
VIBRANT_run.py -t 16 -i vOTU.fasta
./interproscan.sh -i votu.faa -f tsv -dp -goterms -iprlookup -cpu 16 -o my_interproscan/votu_genes.tsv
defense-finder run votu.faa -o defensefinder -a
### host prediction
gtdbtk de_novo_wf --genome_dir /MAG_drep/dereplicated_genomes/ --bacteria --outgroup_taxon p__Patescibacteria --out_dir /host/GTDB-tk_results_1/ --cpus 32 --force --extension fa
gtdbtk de_novo_wf --genome_dir /MAG_drep/dereplicated_genomes/ --archaea --outgroup_taxon p__Altiarchaeota --out_dir /host/GTDB-tk_results_2/ --cpus 32 --force --extension fa
iphop add_to_db --fna_dir /MAG_drep/dereplicated_genomes/ --gtdb_dir /host/GTDB-tk_results/ --out_dir /host/Jun_2025_pub_rw_w_HP_hosts --db_dir /gxfs_work/geomar/smomw651/database/viwrap/ViWrap_db/iPHoP_db/Jun_2025_pub_rw
iphop predict --fa_file /viral/vOTU.fasta --db_dir /hostJun_2025_pub_rw_w_HP_hosts --out_dir /iphop_out -t 32
### diversity
for sample in "$CLEAN_READS_ROOT"/HP_*; do
    sample_id=$(basename "$sample")
    echo "Processing sample: $sample_id"

    # 获取 fastq 文件路径
    r1=$(ls "$sample"/02_CLEAN_READS/*_1.fastq* 2>/dev/null | head -n1)
    r2=$(ls "$sample"/02_CLEAN_READS/*_2.fastq* 2>/dev/null | head -n1)

    if [[ -f "$r1" && -f "$r2" ]]; then
        echo "Using R1: $r1"
        echo "Using R2: $r2"

        bowtie2 -p 16 -x "$INDEX_PATH" -1 "$r1" -2 "$r2" | \
            samtools view -bS - > "$BAM_OUTPUT_DIR/$sample_id.bam"
    else
        echo "ERROR: Missing FastQ files for $sample_id"
    fi
done
metapop --input_samples /virus_metapop/bam/ \
        --reference /virus_metapop/ref_genome/ \
        --norm /virus_metapop/norm.tsv \
        --threads 16

######################
### Ca. H. symbioticus analysis
 for g in `cat genomes.txt`
do
    echo
    echo "Working on $g ..."
    echo
    anvi-gen-contigs-database -f ${g}.fasta \
                              -o HP_${g}.db \
                              --num-threads 4 \
                              -n HP_${g}
done

for g in *.db
do
    anvi-run-hmms -c $g --num-threads 4
    anvi-run-ncbi-cogs -c $g --num-threads 4 --cog-data-dir /database/anvio/cog
    anvi-scan-trnas -c $g --num-threads 4
    anvi-run-pfams -c $g --num-threads 4 --pfam-data-dir /database/anvio
done 
  anvi-display-contigs-stats *db --server-only --port
anvi-profile -c HP_bin.3.db \
             --sample-name HP_bin3 \
             --output-dir HP_bin3 \
             --blank
anvi-gen-genomes-storage -e external-genomes.txt \
                         -o Susy-GENOMES.db

anvi-pan-genome -g Susy-GENOMES.db \
                -n susy \
                -o PAN \
                --num-threads 10               

anvi-compute-genome-similarity --external-genomes external-genomes.txt \
                               --program fastANI \
                               --output-dir ANI \
                               --num-threads 4 \
                               --pan-db Susy/Susy-PAN.db
anvi-merge */PROFILE.db -o Susy-MERGED -c HP_bin.3.db   
anvi-script-add-default-collection -p profiles/Susy-MERGED/PROFILE.db 
  anvi-gen-variability-profile -c HP_bin.3.db \
               -p profiles/Susy-MERGED/PROFILE.db \
               -C DEFAULT \
               -b EVERYTHING \
               --min-occurrence 5 \
               --min-coverage-in-each-sample 20 \
               --include-split-names \
               --quince-mode \
               -o HP_bin3-SNVs.txt
anvi-script-snvs-to-interactive HP_bin3-SNVs.txt -o HP_bin3_snvs 
anvi-get-sequences-for-gene-calls \
    -c Susy-GENOMES.db \
    --get-aa-sequences \
    -o Susy-GENOMES.faa

######################
### GTA tree

iqtree2 -s GTA_TerL_tree_aln_clip.faa \
  -m MFP \
  -B 1000 --bnni \
  -T AUTO \
  -pre TerL_MFP"

##### metatranscriptome
fastp \
  -i Hp27_Mar22_B02_tome_R1.fastq.gz \
  -I Hp27_Mar22_B02_tome_R2.fastq.gz \
  -o Hp27_Mar22_B02_tome_R1.clean.fastq.gz \
  -O Hp27_Mar22_B02_tome_R2.clean.fastq.gz \
  -q 20 \
  -u 30 \
  -l 50 \
  --detect_adapter_for_pe \
  --trim_poly_x \
  -w 8 \
  -h sample.fastp.html \
  -j sample.fastp.json

sortmerna \
  --ref /database/rRNA_databases_v4/smr_v4.3_default_db.fasta \
  --reads HpL1_Jun22_H03_tome_R1.clean.fastq.gz \
  --reads HpL1_Jun22_H03_tome_R2.clean.fastq.gz \
  --paired_in \
  --fastx \
  --other HpL1_Jun22_H03.non_rRNA \
  --aligned HpL1_Jun22_H03.rRNA \
  --threads 8 \
  --workdir sortmerna_work/HpL1_Jun22_H03

bowtie2 \
  -x "$project/host_remove/host_index" \
  -1 "$r1" \
  -2 "$r2" \
  --un-conc-gz "$project/host_removed/${base}_nohost_R%.fq.gz

salmon quant \
  -i /salmon_index \
  -l A \
  -1 /SRR34071591_nohost_R1.fq.gz \
  -2 /SRR34071591_nohost_R2.fq.gz \
  -p 8 \
  --meta \
  -o /salmon/SRR34071591

