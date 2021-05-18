import os

configfile: "config.yaml"

data_dir = 'data'
library_dir = 'play_data_ref_annot'
reads = ['1', '2']

fastqc_dir = 'fastqc'
multiqc_dir = 'multiqc'
trimmed_dir = 'trimmed'
genome_dir = 'genome'
star_dir = 'star_aligned'
samtools_dir = 'samtools_sorted'
fcounts_dir = 'fcounts'
deseq_dir = 'DE_results'


rule all:
    input:
        #expand(os.path.join(fastqc_dir, '{sample}_R{read}_001_fastqc.html'), sample=config['samples'], read=reads),
        os.path.join(multiqc_dir, 'multiqc_report.html'),
        #expand(os.path.join(trimmed_dir, '{sample}_R1_001.fastq'), sample=config['samples']),
        #expand(os.path.join(trimmed_dir, '{sample}_R2_001.fastq'), sample=config['samples']),
        #expand(os.path.join(fastqc_dir, trimmed_dir, '{sample}_R{read}_001_fastqc.html'), sample=config['samples'], read=reads),
        os.path.join(multiqc_dir, trimmed_dir, 'multiqc_report.html'),
        #os.path.join(genome_dir, 'genomeParameters.txt'),
        #os.path.join(genome_dir, 'Genome'),
        #expand(os.path.join(star_dir, '{sample}', 'Aligned.sortedByCoord.out.bam'), sample=config['samples']),
        #expand(os.path.join(samtools_dir, '{sample}.bam'), sample=config['samples']),
        os.path.join(fcounts_dir, 'fcount_result.txt'),
        expand(os.path.join(deseq_dir, '{prep}_DE_volcano.pdf'), prep=config['preps']),
        expand(os.path.join(deseq_dir, '{prep}_DE_results.csv'), prep=config['preps'])
               
                
rule fastqc:
    input:
       os.path.join(data_dir, '{sample}_R{read}_001.fastq')
    output:
        os.path.join(fastqc_dir, '{sample}_R{read}_001_fastqc.html')
    shell:
        "fastqc {input} --outdir={fastqc_dir}"


rule multiqc: 
    input: 
        expand(os.path.join(fastqc_dir, '{sample}_R{read}_001_fastqc.html'), sample=config['samples'], read=reads)
    output:
        os.path.join(multiqc_dir, 'multiqc_report.html')
    shell: 
        "multiqc {fastqc_dir} -o {multiqc_dir}"
        
        
rule bbduk:
    input:
        r1 = os.path.join(data_dir, '{sample}_R1_001.fastq'),
        r2 = os.path.join(data_dir, '{sample}_R2_001.fastq')
    output:
        out1 = os.path.join(trimmed_dir, '{sample}_R1_001.fastq'),
        out2 = os.path.join(trimmed_dir, '{sample}_R2_001.fastq')
    shell: 
        "scripts/bbmap/bbduk.sh in1={input.r1} out1={output.out1} in2={input.r2} out2={output.out2} ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10"
        
        
rule fastqc_trimmed:
    input:
       os.path.join(trimmed_dir, '{sample}_R{read}_001.fastq')
    output:
        os.path.join(fastqc_dir, trimmed_dir, '{sample}_R{read}_001_fastqc.html')
    shell:
        "fastqc {input} --outdir={fastqc_dir}/{trimmed_dir}"


rule multiqc_trimmed: 
    input: 
        expand(os.path.join(fastqc_dir, trimmed_dir, '{sample}_R{read}_001_fastqc.html'), sample=config['samples'], read=reads)
    output:
        os.path.join(multiqc_dir, trimmed_dir, 'multiqc_report.html')
    shell: 
        "multiqc {fastqc_dir}/{trimmed_dir} -o {multiqc_dir}/{trimmed_dir}"
        
        
rule star_generate:
    input:
        gtf = os.path.join(library_dir, 'chr19_20Mb.gtf'),
        fa = os.path.join(library_dir, 'chr19_20Mb.fa')
    output:
        os.path.join(genome_dir, 'genomeParameters.txt'),
        os.path.join(genome_dir, 'Genome')
    shell:
        "scripts/STAR --genomeDir {genome_dir} --runMode genomeGenerate --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --genomeSAindexNbases 11 --outFileNamePrefix star_files"
        
        
rule star_align:
    input:
        os.path.join(genome_dir, 'genomeParameters.txt'),
        os.path.join(genome_dir, 'Genome'),
        r1 = os.path.join(trimmed_dir, '{sample}_R1_001.fastq'),
        r2 = os.path.join(trimmed_dir, '{sample}_R2_001.fastq'),
        gtf = os.path.join(library_dir, 'chr19_20Mb.gtf')       
    output:  
        os.path.join(star_dir, '{sample}', 'Aligned.sortedByCoord.out.bam')  
    shell:    
        "scripts/STAR --readFilesIn {input.r1} {input.r2} --genomeDir {genome_dir} --outSAMtype BAM SortedByCoordinate --sjdbGTFfile {input.gtf} --outFileNamePrefix {star_dir}/{wildcards.sample}/"
        
        
rule samtools_sort:
    input:
        os.path.join(star_dir, '{sample}', 'Aligned.sortedByCoord.out.bam')
    output:
        os.path.join(samtools_dir, '{sample}.bam')
    shell:
        "samtools sort -n -T sorted_reads/{wildcards.sample} -O bam {input} > {output}"


def fcounts_inputs(wildcards):
    files = expand(os.path.join(samtools_dir, '{sample}.bam'), sample=config['samples'])
    return files


rule feature_counts:
    input:
        fcounts_inputs
    output:
        os.path.join(fcounts_dir, 'fcount_result.txt')
    shell:
        "scripts/featureCounts {input} -p -t exon -g gene_id -a {library_dir}/chr19_20Mb.gtf -o {output} -s 1"
        
        
rule de_analysis:
    input:
        os.path.join(fcounts_dir, 'fcount_result.txt'),
    params:
        deseq_dir,
        config['preps']
    output:
        os.path.join(deseq_dir, '{prep}_DE_volcano.pdf'),
        os.path.join(deseq_dir, '{prep}_DE_results.csv')
    script:
        "scripts/deseq2.R"
