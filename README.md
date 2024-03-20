# LRAPmut
Weeds is a germline small variant caller for hybrid of long and short reads. 

# Install
```
conda create -n Weeds
conda install -c bioconda bcftools samtools clair3=1.0.4 python=3.9.0 -y
git clone https://github.com/Roick-Leo/Weeds.git
chmod -R a+x ./Weeds
```

# Usage
## General Usage
```
./lrapmut.bin \
    --env_dir ${The_bin_dir_of_Weeds_conda_env} \
    --dnb_bam_fn ${DNB_BAM} \
    --cyclone_bam_fn ${CYCLONE_bam} \
    --ref_fn ${REF} \
    --threads ${THREADS} \ 
    --model_path ${MODEL_DIR} \
    --output ${OUTPUT_DIR} \
    --longreads_region ${LONGREADS_BED}
```

## Demo
```
minimap2 -ax map-ont --secondary=no -t 15 /path/to/ref /path/to/fastq | samtools view -@ 15 -bS | samtools sort -@ 15 -o /path/to/output_sorted.bam

/path/to/Weeds.bin --env_dir /path/to/Weeds_conda_env_dir/bin --dnb_bam_fn /path/to/dnb_bam_file --cyclone_bam_fn /path/to/cyclone_bam_file --ref_fn /path/to/GRch38_ref.fa --model_path /path/to/model_dir 
--longreads_region /path/to/Weeds/data/Longreads.bed --threads 15 --output /path/to/output_dir
```

# Final comparison
```
mkdir -p /path/to/hppy_output_dir
hap.py \
/path/to/benchmark/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
/path/to/LRAPmut/Final_merged.vcf.gz \
-f /path/to/benchmark/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
-r /path/to/reference/hg38_noalt_withrandom/hg38.fa \
-o /path/to/hppy_output_dir \
--engine=vcfeval \
--threads=15 \
--pass-only
```