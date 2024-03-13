import argparse
import os
import subprocess
from datetime import datetime

parser = argparse.ArgumentParser(description="Mutation screening process combining nanopore sequencing and next-generation sequencing.")
parser.add_argument('--env_dir', required=True, type=str,
                    help="Path to conda env.")
parser.add_argument('--bed_fn', type=str, default="",
                    help="Call variants only in the provided bed regions.")
parser.add_argument('--dnb_bam_fn', required=True, type=str,
                    help="BAM file input of DNB. The input file must be samtools indexed.")
parser.add_argument('--cyclone_bam_fn', required=True, type=str,
                    help="BAM file input of cyclone. The input file must be samtools indexed.")
parser.add_argument('--ref_fn', required=True, type=str,
                    help="FASTA reference file input. The input file must be samtools indexed.")
parser.add_argument('--model_path', required=True, type=str,
                    help="The folder path containing the model (requiring six files in the folder, including pileup.data-00000-of-00002, pileup.data-00001-of-00002 pileup.index, full_alignment.data-00000-of-00002, full_alignment.data-00001-of-00002 and full_alignment.index).")
parser.add_argument("--longreads_region", required=True,
                    help="The f1-score was confirmed to be higher in these regions.")
parser.add_argument("--out_dir", required=True,
                    help="The output dir.")
parser.add_argument("--sample", required=True,
                    help="Output prefix.")
parser.add_argument("--bcftools", default="bcftools", type=str,
                    help="The path of bcftools.")
parser.add_argument("--tabix", default="tabix", type=str,
                    help="The path of tabix.")
parser.add_argument("--threads", default=4, type=int,
                    help="Threads to be used.")
args = parser.parse_args()

def vcf_iter(vcf_file):
    with open(vcf_file,"r") as vcf:
        for line in vcf:
            if not line.startswith("#"):
                yield line
chrom_to_num = dict(zip(["chr"+str(i) for i in range(1,23)]+["chrX","chrY"],list(range(1,25))))
def get_rec_list(rec):
    ll = rec.strip().split("\t")
    QUAL = float(ll[5])
    CHROM = ll[0]
    REF = ll[3]
    if CHROM not in ["chrX","chrY"]:
        CHROM = int(CHROM.replace("chr",""))
    elif CHROM == "chrX":
        CHROM = 23
    else:
        CHROM = 24
    AF = ll[-1].split(":")[-1]
    ALT = ll[4]
    if "," in AF:
        AF_ll = [float(i) for i in AF.split(",")]
        ALT_ll = ALT.split(",")
        AF_tmp = AF_ll[0]
        ALT_tmp = ALT_ll[0]
        for i in range(1,len(AF_ll)):
            if AF_ll[i] > AF_tmp:
                AF_tmp = AF_ll[i] 
                ALT_tmp = ALT_ll[i]
        AF = AF_tmp
        ALT = ALT_tmp
    else:
        AF = float(AF)
    POS = int(ll[1])
    return [CHROM,POS,REF,ALT,QUAL,AF]
def snp_calling(TOOL,bam,ref,threads,platform,model,output,sample,bed):
    cmd="{} --bam_fn {} --ref_fn {} --threads {} --platform {} --model_path {} --output{} --sample_name {}"
    if bed != "":
        cmd += " --bed_fn {}".format(bed)
    subprocess.run(cmd.format(TOOL,bam,ref,threads,platform,model,output,sample),shell=True, stdout=subprocess.DEVNULL)

conda_env = args.env_dir
TOOL = os.path.join(conda_env,"run_clair3.sh")
DNB_bam = args.dnb_bam_fn
Cyclone_bam = args.cyclone_bam_fn
Model_path = args.model_path
DNB_model = os.path.join(Model_path,"DNB")
Cyclone_model = os.path.join(Model_path,"Cyclone")
REF = args.ref_fn
BED = args.bed_fn

# Make output dir
out_dir = args.out_dir
DNB_tmp_dir = os.path.join(out_dir,"DNB")
Cyclone_tmp_dir = os.path.join(out_dir,"Cyclone")
if not os.path.exists(out_dir):
    subprocess.run("mkdir -p {} {} {}".format(out_dir,DNB_tmp_dir,Cyclone_tmp_dir))

c_vcf = os.path.join(Cyclone_tmp_dir,"merge_output.vcf.gz")
d_vcf = os.path.join(DNB_tmp_dir,"merge_output.vcf.gz")
BCFTOOLS = args.bcftools
threads = args.threads
lone_reads_bed_file = args.longreads_region
sample = args.sample

## SNP calling
# Cyclone
now = datetime.now()
formatted_time = now.strftime("%Y-%m-%d %H:%M:%S")
print("[{}] STEP1 Call variants using cyclone data".format(formatted_time))
snp_calling(TOOL,Cyclone_bam,REF,threads,"ont",Cyclone_model,Cyclone_tmp_dir,sample,BED)
# DNB
now = datetime.now()
formatted_time = now.strftime("%Y-%m-%d %H:%M:%S")
print("[{}] STEP2 Call variants using DNB data".format(formatted_time))
snp_calling(TOOL,DNB_bam,REF,threads,"ont",DNB_model,DNB_tmp_dir,sample,BED)

## Get unzip vcf
# Cyclone
if c_vcf.endswith(".vcf.gz"):
    c_new_vcf = c_vcf.replace(".vcf.gz",".vcf")
    subprocess.run("{} view --threads {} -Ov -o {} {}".format(BCFTOOLS,threads,c_new_vcf,c_vcf),shell=True, stdout=subprocess.DEVNULL)
    c_vcf = c_new_vcf
else:
    pass
# DNB
if d_vcf.endswith(".vcf.gz"):
    d_new_vcf = d_vcf.replace(".vcf.gz",".vcf")
    subprocess.run("{} view --threads {} -Ov -o {} {}".format(BCFTOOLS,threads,d_new_vcf,d_vcf),shell=True, stdout=subprocess.DEVNULL)
    d_vcf = d_new_vcf
else:
    pass

dnb_iter = vcf_iter(d_vcf)
nano_iter = vcf_iter(c_vcf)
bed_iter = iter(open(lone_reads_bed_file,"r"))
dnb_rec = next(dnb_iter,None)
nano_rec = next(nano_iter,None)
region = next(bed_iter,None)
output = os.path.join(out_dir,sample+".vcf")
with open(c_vcf,"r") as cv, open(output, "w") as out:
    for line in cv:
        if line.startswith("#"):
            out.write(line)
        else:
            break
    out.close()

with open(output,"a") as out_vcf:
    while dnb_rec is not None and nano_rec is not None:
        dnb_chrom,dnb_pos,dnb_ref,dnb_alt,dnb_qual,dnb_af = get_rec_list(dnb_rec)
        if len(dnb_alt) != len(dnb_ref):
            out_vcf.write(str(dnb_rec))
            dnb_rec = next(dnb_iter, None)
            continue

        nano_chrom,nano_pos,nano_ref,nano_alt,nano_qual,nano_af = get_rec_list(dnb_rec)
        if nano_af < 0.1:
            nano_rec = next(nano_iter, None)
            continue
        if len(nano_alt) != len(nano_ref):
            nano_rec = next(nano_iter, None)
            continue
        
        if region is not None:
            region_list = region.strip().split("\t")
            region_chrom = region_list[0]
            region_chrom = chrom_to_num[region_chrom]
            region_start = int(region_list[1])
            region_end = int(region_list[2])
        else:
            region_chrom = 24
            region_start = 10**100
            region_end = 10**100
        
        ## merge
        if dnb_chrom == nano_chrom:
            if dnb_pos == nano_pos:
                if dnb_alt == nano_alt:
                    out_vcf.write(str(dnb_rec))
                elif (dnb_chrom > region_chrom) or (dnb_chrom == region_chrom and dnb_pos > region_end):
                    region = next(bed_iter,None)
                    continue
                elif ((dnb_chrom < region_chrom) or (dnb_chrom == region_chrom and dnb_pos < region_start)):
                    if nano_qual > dnb_qual and (nano_qual >= 22 or (nano_qual >= 10 and nano_af >= 0.8)):
                        out_vcf.write(str(nano_rec))
                    else:
                        out_vcf.write(str(dnb_rec))
                else:
                    if dnb_qual > nano_qual and (dnb_qual >= 20 or (dnb_qual >= 10 and dnb_af >= 0.8)):
                        out_vcf.write(str(dnb_rec))
                    else:
                        out_vcf.write(str(nano_rec))
                dnb_rec = next(dnb_iter,None)
                nano_rec = next(nano_iter,None)
                continue
            elif dnb_pos < nano_pos:
                if (dnb_chrom > region_chrom) or (dnb_chrom == region_chrom and dnb_pos > region_end):
                    region = next(bed_iter,None)
                    continue
                elif ((dnb_chrom < region_chrom) or (dnb_chrom == region_chrom and dnb_pos < region_start)):
                    out_vcf.write(str(dnb_rec))
                elif dnb_qual >= 20 or (dnb_qual >= 10 and dnb_af >= 0.8):
                    out_vcf.write(str(dnb_rec))
                dnb_rec = next(dnb_iter,None)
                continue
            else:
                if (nano_chrom > nano_chrom) or (nano_chrom == region_chrom and nano_pos > region_end):
                    region = next(bed_iter,None)
                    continue
                elif ((nano_chrom < region_chrom) or (nano_chrom == region_chrom and nano_pos < region_start)):
                    if nano_qual >= 22 or (nano_qual >= 10 and nano_af >= 0.8):
                        out_vcf.write(str(nano_rec))
                else:
                    out_vcf.write(str(nano_rec))
                nano_rec = next(nano_iter,None)
                continue
        elif dnb_chrom < nano_chrom:
            if (dnb_chrom > region_chrom) or (dnb_chrom == region_chrom and dnb_pos > region_end):
                region = next(bed_iter,None)
                continue
            elif ((dnb_chrom < region_chrom) or (dnb_chrom == region_chrom and dnb_pos < region_start)):
                out_vcf.write(str(dnb_rec))
            elif dnb_qual >= 20 or (dnb_qual >= 10 and dnb_af >= 0.8):
                out_vcf.write(str(dnb_rec))
            dnb_rec = next(dnb_iter,None)
            continue
        else:
            if (nano_chrom > nano_chrom) or (nano_chrom == region_chrom and nano_pos > region_end):
                region = next(bed_iter,None)
                continue
            elif ((nano_chrom < region_chrom) or (nano_chrom == region_chrom and nano_pos < region_start)):
                if nano_qual >= 22 or (nano_qual >= 10 and nano_af >= 0.8):
                    out_vcf.write(str(nano_rec))
            else:
                out_vcf.write(str(nano_rec))
            nano_rec = next(nano_iter,None)
            continue

    while dnb_rec is not None:
        ddnb_chrom,dnb_pos,dnb_ref,dnb_alt,dnb_qual,dnb_af = get_rec_list(dnb_rec)
        if len(dnb_alt) != len(dnb_ref):
            out_vcf.write(str(dnb_rec))
            dnb_rec = next(dnb_iter, None)
            continue
        out_vcf.write(str(dnb_rec))
        dnb_rec = next(dnb_iter,None)

    while nano_rec is not None:
        nano_chrom,nano_pos,nano_ref,nano_alt,nano_qual,nano_af = get_rec_list(dnb_rec)
        if nano_af < 0.1 or nano_qual < 1:
            nano_rec = next(nano_iter, None)
            continue
        if len(nano_alt) != len(nano_ref):
            nano_rec = next(nano_iter, None)
            continue
        if nano_qual >= 22 or (nano_qual >= 10 and nano_af >= 0.8):
            out_vcf.write(str(nano_rec))
        nano_rec = next(nano_iter, None)

TABIX = args.tabix
subprocess.run("{} view --threads {} {} -Oz -o {}; {} {} ; rm {} {}".format(BCFTOOLS,threads,output,output+".gz",TABIX,output+".gz",c_vcf,d_vcf),shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)