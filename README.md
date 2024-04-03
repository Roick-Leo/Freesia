# Freesia
Freesia is a germline small variant caller for hybrid of long and short reads. 

# Install
## Step1. Configure the conda source
```
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
conda config --set show_channel_urls yes
```
## Step2. Install Freesia
```
conda create -n Freesia
conda activate Freesia
conda install -c bioconda -c conda-forge bcftools=1.17 samtools clair3=1.0.4 python=3.9.0 -y
git clone https://github.com/Roick-Leo/Freesia.git
chmod -R a+x ./Freesia
```
## Step3. Install Boost Graph Library for realignment process
```
conda activate Freesia
conda install -c conda-forge boost -y
# CONDA_PREFIX can be obtained by running conda env list
CONDA_PREFIX="/path/to/conda/envs/Freesia"

# Optional steps: update gcc & g++
conda install -c conda-forge gcc libgcc gxx_linux-64 -y
cd ${CONDA_PREFIX}/bin
ln -s x86_64-conda_cos6-linux-gnu-g++ g++ 

cd ${CONDA_PREFIX}/bin/preprocess/realign
g++ -std=c++14 -O1 -shared -fPIC -o realigner ssw_cpp.cpp ssw.c realigner.cpp
g++ -std=c++11 -shared -fPIC -o debruijn_graph -O3 debruijn_graph.cpp -I ${CONDA_PREFIX}/include -L ${CONDA_PREFIX}/lib
```

# Usage
## General Usage
```
./lrapmut.bin \
    --env_dir ${The_bin_dir_of_Freesia_conda_env} \     # the bin dir of freesia environment
    --dnb_bam_fn ${DNB_BAM} \                           # the sorted bam file of DNB-seq with index
    --cyclone_bam_fn ${CYCLONE_bam} \                   # the sorted bam file of Cyclone with index
    --ref_fn ${REF} \                                   # the path to refrence file
    --threads ${THREADS} \                              # thr threads num
    --model_path ${MODEL_DIR} \                         # the path to model
    --output ${OUTPUT_DIR} \                            # the path to output dir
    --longreads_region ${LONGREADS_BED}                 # the path to bed file where is available in the project document
```

## Demo
```
# minimap2
minimap2 -ax map-ont --secondary=no -t 15 /path/to/ref /path/to/fastq | samtools view -@ 15 -bS | samtools sort -@ 15 -o /path/to/output_sorted.bam

# bwa
bwa mem -t 15 -m 50G -R '@RG\tID:foo_lane\tPL:BGI\tLB:library\tSM:HG002_pcr_free' /path/to/ref /path/to/fastq_1 /path/to/fastq_2 -M -Y -o /path/to/output.sam
samtools view -bS /path/to/output.sam | samtools sort -@ 15 -o /path/to/output_sorted.bam
samtools index /path/to/output_sorted.bam

/path/to/Freesia.bin --env_dir /path/to/Freesia_conda_env_dir/bin --dnb_bam_fn /path/to/dnb_bam_file --cyclone_bam_fn /path/to/cyclone_bam_file --ref_fn /path/to/GRch38_ref.fa --model_path /path/to/Freesia/model_dir 
--longreads_region /path/to/Freesia/data/Longreads.bed --threads 15 --output /path/to/output_dir
```

# Final comparison
```
mkdir -p /path/to/hppy_output_dir
hap.py \
/path/to/benchmark/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
/path/to/Freesia/Freesia_out.vcf.gz \
-f /path/to/benchmark/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
-r /path/to/reference/hg38_noalt_withrandom/hg38.fa \
-o /path/to/hppy_output_dir \
--engine=vcfeval \
--threads=15 \
--pass-only
```

# Contact
- He Lei
  - Email: helei1@genomics.cn