import re
import os

configfile: "config.yaml"


DIR_OUT = os.path.join(config["path"]["dir_out"],"")
DIR_BAM = os.path.join(config["path"]["dir_bam"],"")
DIR_STATS = os.path.join(config["path"]["dir_stats"],"")
DIR_LOG = os.path.join(config["path"]["dir_log"],"")

#(wholenames,) = glob_wildcards(DIR_FASTQ+"{wholename}.fastq.gz")
(wholenames,) = glob_wildcards(DIR_BAM+"{wholename}.bam")

profiletypes = config["summary"]["profiletypes"]
BINSIZES=config["QDNAseq"]["BINSIZES"]
imagetype=config["ACE"]["imagetype"]
ACEBINSIZES=config["ACE"]["ACEBINSIZES"]
setting = config["pipeline"]["setting"]
REF_CLONALITY = config["Clonality"]["reference"]

def getnames():
    SAMPLES=dict()
    for wholename in wholenames:
        sample = wholename
        bamfile=DIR_BAM+wholename+".bam"
        SAMPLES[sample]=bamfile
    return(SAMPLES)

SAMPLES=getnames()
print(SAMPLES)

#============================================Run_setting====================================================
if setting == "service": #rule service
    rule service:
        input:
            expand(DIR_OUT + "{binSize}kbp/summary.html", binSize=BINSIZES), # summary
elif setting == "test_run": # FOR A TEST RUN, STOPS AT rule QDNAseq_segment !!!
    rule test_run:
        input:
            expand(DIR_OUT + "{binSize}kbp/data/Clonality/{binSize}.report.png",binSize=BINSIZES),
            #expand(DIR_OUT + "{binSize}kbp/ACE/{ploidy}N/segmentfiles/{sample}_segments.tsv", binSize=ACEBINSIZES, ploidy=config["ACE"]["ploidies"], sample=SAMPLES.keys()),
            #expand(DIR_OUT + "{binSize}kbp/data/{binSize}kbp-called.rds",binSize=BINSIZES)
elif setting == "research": #rule research
    rule research:    
        input:
            expand(DIR_OUT + "{binSize}kbp/summary.html", binSize=BINSIZES), # summary
            #expand(DIR_OUT + "{binSize}kbp/ACE/{ploidy}N/segmentfiles/{sample}_segments.tsv", binSize=ACEBINSIZES, ploidy=config["ACE"]["ploidies"], sample=SAMPLES.keys()), # postanalysisloop_ACE
            #expand(DIR_OUT + "{ACEbinSize}kbp/profiles/call_cellularity_based/index.html",ACEbinSize=ACEBINSIZES), # lightBox
else: #rule all
    rule all:
        input: 
            expand(DIR_OUT + "{binSize}kbp/summary.html", binSize=BINSIZES), # summary
            #expand(DIR_OUT + "{binSize}kbp/ACE/{ploidy}N/segmentfiles/{sample}_segments.tsv", binSize=ACEBINSIZES, ploidy=config["ACE"]["ploidies"], sample=SAMPLES.keys()), # postanalysisloop_ACE
            #expand(DIR_OUT + "{ACEbinSize}kbp/data/{ACEbinSize}kbp-call_cellularity_based.rds", ACEbinSize=ACEBINSIZES), # CNA_call_cellularity_based
            expand(DIR_OUT + "{ACEbinSize}kbp/profiles/call_cellularity_based/index.html",ACEbinSize=ACEBINSIZES), # lightBox
            ##expand(DIR_OUT + "{binSize}kbp/ACE/{ploidy}N/{sample}/summary_{sample}.{imagetype}", imagetype=imagetype ,binSize=ACEBINSIZES, ploidy=config["ACE"]["ploidies"], sample=SAMPLES.keys()),

#============================================Rules=======================================================

rule generate_stats:
    input:
        bam=DIR_OUT + DIR_BAM + "{sample}.bam",
        script="scripts/Run_generate_stats.sh"
    output:
        DIR_OUT + DIR_STATS + "{sample}.reads.all"
    params:
        outdir=DIR_OUT + DIR_STATS
    shell:
        "{input.script} {input.bam} {wildcards.sample} {params.outdir} {output}"

rule QDNAseq_binReadCounts:
    input:
        bams=expand(DIR_BAM + "{sample}.bam", sample=SAMPLES.keys()),
    output:
        binReadCounts=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-raw.rds"
    params:
        genome=config["QDNAseq"]["genome"],
        suppressMessages=config["pipeline"]["suppressMessages"]
    log: DIR_OUT + DIR_LOG + "QDNAseq/{binSize}kbp/binReadCounts.log"
    script:
        "scripts/Run_QDNAseq_binReadCounts.R"

rule QDNAseq_normalize:
    input:
        binReadCounts=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-raw.rds",
    output:
        corrected=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-corrected.rds",
        allprofiles=expand(DIR_OUT + "{{binSize}}kbp/profiles/corrected/{samples}.png",samples=SAMPLES.keys()),
    params:
        profiles=DIR_OUT + "{binSize}kbp/profiles/corrected/",
        chrom_filter=config["QDNAseq"]["chrom_filter"],
        suppressMessages=config["pipeline"]["suppressMessages"]
    log: DIR_OUT + DIR_LOG + "QDNAseq/{binSize}kbp/normalizeBins.log"
    script:
        "scripts/Run_QDNAseq_normalize.R"

rule deWave:
    input:
        corrected=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-corrected.rds",
    output:
        dewaved=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-dewaved.rds",
        allprofiles=expand(DIR_OUT + "{{binSize}}kbp/profiles/dewaved/{samples}.png",samples=SAMPLES.keys()),
    params:
        profiles=DIR_OUT + "{binSize}kbp/profiles/dewaved/",
        dewave_data=config["QDNAseq"]["dewave_data"],
        suppressMessages=config["pipeline"]["suppressMessages"]
    log: DIR_OUT + DIR_LOG + "QDNAseq/{binSize}kbp/dewave.log"
    script:
        "scripts/Run_deWave.R"

rule QDNAseq_segment:
    input:
        dewaved=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-dewaved.rds" if config['QDNAseq']['dewave'] else DIR_OUT + "{binSize}kbp/data/{binSize}kbp-corrected.rds",
    output:
        segmented=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-segmented.rds",
        allprofiles=expand(DIR_OUT + "{{binSize}}kbp/profiles/segmented/{samples}.png",samples=SAMPLES.keys()),
        copynumbers=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-copynumbers.igv",
        segments=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-segments.igv",
        copynumbersbed=expand(DIR_OUT + "{{binSize}}kbp/BED/{samples}-copynumbers.bed",samples=SAMPLES.keys()),
        segmentsbed=expand(DIR_OUT + "{{binSize}}kbp/BED/{samples}-segments.bed",samples=SAMPLES.keys()),
    params:
        profiles=DIR_OUT + "{binSize}kbp/profiles/segmented/",
        failed=DIR_OUT + "{binSize}kbp/failed_samples.txt",
        minimal_used_reads=config["QDNAseq"]["minimal_used_reads"],
        copynumbersbed=DIR_OUT + "{binSize}kbp/BED/%s-copynumbers.bed",
        segmentsbed=DIR_OUT + "{binSize}kbp/BED/%s-segments.bed",
        bedfolder=DIR_OUT + "{binSize}kbp/BED/",
        suppressMessages=config["pipeline"]["suppressMessages"]
    log: DIR_OUT + DIR_LOG + "QDNAseq/{binSize}kbp/segment.log"
    script:
        "scripts/Run_QDNAseq_segment.R"

#----------------------------------------------------------FOR_LUNG_CANCER----------------------------------------------------------
rule clonality:
    input:
        RDS=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-segmented.rds"
    output:
        profiles=directory(DIR_OUT + "{binSize}kbp/data/Clonality/"),
        report_png=DIR_OUT + "{binSize}kbp/data/Clonality/{binSize}.report.png",
        report_bmp=DIR_OUT + "{binSize}kbp/data/Clonality/{binSize}.report.bmp"
    params:
        suppressMessages=config["pipeline"]["suppressMessages"],
        projname="{binSize}-Clonality",
        reference=REF_CLONALITY
    log: DIR_OUT + DIR_LOG + "QDNAseq/{binSize}kbp/clonality.log"
    script:
        "scripts/d4-pipeline-master/bin/clonality.R"

#----------------------------------------------------------FOR_MELANOMA----------------------------------------------------------
rule ACE:
    input:
        segmented=DIR_OUT + "{ACEbinSize}kbp/data/{ACEbinSize}kbp-segmented.rds",
    output:
        #ACE=expand(DIR_OUT + "{{ACEbinSize}}kbp/ACE/{{ploidy}}N/summary_files/summary_{sample}.{imagetype}", sample=SAMPLES.keys()),
        fitpicker=DIR_OUT + "{ACEbinSize}kbp/ACE/{ploidy}N/fitpicker_{ploidy}N.tsv",
    params:
        outputdir=DIR_OUT + "{ACEbinSize}kbp/ACE/",
        failed=DIR_OUT + "{ACEbinSize}kbp/failed_samples.txt",
        suppressMessages=config["pipeline"]["suppressMessages"]
    log:DIR_OUT + DIR_LOG + "ACE/{ACEbinSize}kbp/{ploidy}N/ACE_log.tsv"
    script:
        "scripts/Run_ACE.R"

rule postanalysisloop_ACE:
    input:
        segmented=DIR_OUT + "{ACEbinSize}kbp/data/{ACEbinSize}kbp-segmented.rds",
        fitpicker=DIR_OUT + "{ACEbinSize}kbp/ACE/{ploidy}N/fitpicker_{ploidy}N.tsv",
    output:
        ACE_post=expand(DIR_OUT + "{{ACEbinSize}}kbp/ACE/{{ploidy}}N/segmentfiles/{sample}_segments.tsv", sample=SAMPLES.keys()),
    params:
        outputdir=DIR_OUT + "{ACEbinSize}kbp/ACE/{ploidy}N/",
        failed=DIR_OUT + "{ACEbinSize}kbp/failed_samples.txt",
        suppressMessages=config["pipeline"]["suppressMessages"]
    log:DIR_OUT + DIR_LOG + "ACE/{ACEbinSize}kbp/{ploidy}N/ACE_post_log.tsv"
    script:
        "scripts/Run_postanalysisloop_ACE.R"

rule CNA_call:
    input:
        segmented=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-segmented.rds",
    output:
        called=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-called.rds",
        freqplot=DIR_OUT + "{binSize}kbp/profiles/freqPlot/freqPlot_{binSize}kbp.png",
        allprofiles=expand(DIR_OUT + "{{binSize}}kbp/profiles/called/{samples}.png",samples=SAMPLES.keys()),
        calls=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-calls.igv"
    params:
        profiles=DIR_OUT + "{binSize}kbp/profiles/called/",
        failed=DIR_OUT + "{binSize}kbp/failed_samples.txt",
        suppressMessages=config["pipeline"]["suppressMessages"]
    log: DIR_OUT + DIR_LOG + "CNA/{binSize}kbp/call.log"
    script:
        "scripts/Run_CNA_call.R"

'''
rule lightBox:
    input:
        sample=expand(DIR_OUT + "{{binSize}}kbp/profiles/{{profiletype}}/{sample}.png", sample=SAMPLES.keys()),
        script="scripts/createLightBox.sh",
    output:
        index=DIR_OUT + "{binSize}kbp/profiles/{profiletype}/index.html",
    params:
        profiles=DIR_OUT + "{binSize}kbp/profiles/{profiletype}/",
        lb2dir="lb2/"
    shell:
        "{input.script} {params.profiles} {params.lb2dir} > {output.index}"

rule summary:
#TODO ADAPT to start from bam files.
    input:
        stats=expand(DIR_OUT + DIR_STATS + "{sample}.reads.all", sample=SAMPLES.keys()),
        index=expand(DIR_OUT + "{{binSize}}kbp/profiles/{profiletype}/index.html", profiletype=profiletypes),
        script="scripts/Run_lane-summary.sh",
        #qcfastq=expand(DIR_OUT + DIR_QC +  "qc-fastq/{wholename}_fastqc.html", wholename=wholenames),
        #bamqc=expand(DIR_OUT + DIR_QC +  "qc-bam/{sample}_fastqc.html", sample=SAMPLES.keys()),
    output:
        DIR_OUT + "{binSize}kbp/summary.html"
    params:
        bamfolder=DIR_OUT + DIR_BAM,
        statsfolder = DIR_OUT + DIR_STATS,
        #qcFastqfolder= DIR_OUT + DIR_QC + "qc-fastq/",
        #qcBamfolder= DIR_OUT + DIR_QC + "qc-bam/"
    shell:
        "{input.script} {wildcards.binSize}kpb {params.bamfolder} {params.statsfolder} {params.qcFastqfolder} {params.qcBamfolder}> {output}"

'''
#----------------------------------------------------------------------------------------------------------------
