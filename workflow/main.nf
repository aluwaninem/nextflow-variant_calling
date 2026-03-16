#!/usr/bin/env nextflow

workflow {
   def fastq = Channel.fromFilePairs(params.fastq)
   def genome = Channel.fromPath(params.chr19Folder)

   FASTQC(fastq)
   TRIMMOMATIC(fastq)
   BWA(genome,TRIMMOMATIC.out.trimmed)
   SAM_TO_BAM(BWA.out.bwa_sam)
   VAR_CAL(SAM_TO_BAM.out.bam_file, genome)

}

process FASTQC {   
  publishDir "${params.outdir}/fastqc", mode: 'copy', overwrite:true

   input:
   tuple val(id), path(reads)

   output:
   tuple val(id), path("*_fastqc.*")
   
   script:
   """
   module load fastqc
   fastqc $reads
   """
}

process TRIMMOMATIC {
  publishDir "${params.outdir}/trimmomatic", mode: 'copy', overwrite:true

   input:
   tuple val(id), path(reads)

   output:
   tuple val(id), path("${id}_1.trimmed.fastq"), path("${id}_2.trimmed.fastq"), emit: trimmed

   script:
   """
   singularity exec ${params.sif} java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
   -threads ${task.cpus} -phred33 \
   ${reads[0]} ${reads[1]} \
   ${id}_1.trimmed.fastq ${id}_1.unpaired.fastq \
   ${id}_2.trimmed.fastq ${id}_2.unpaired.fastq \
   LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:32
   """ 
}


process BWA {

   input:
   path chr19
   tuple val(id), path(read1), path(read2)

   output:
   tuple val(id), path("*.sam"), emit:bwa_sam

   script:
   """
   module load bwa
   bwa mem -t 2 ${chr19}/chr19.fa $read1 $read2 > ${id}.sam
   """
}

process SAM_TO_BAM { 

   input:
   tuple val(id), path(sam)

   output:
   tuple val(id), path("*-sorted.bam"), path("*-sorted.bam.bai"), emit:bam_file

   script:
   """
   module load samtools
   
   samtools view -b $sam -o ${id}.bam
   samtools sort -m 5G -@ 2 -T tmp -o ${id}-sorted.bam ${id}.bam
   samtools index ${id}-sorted.bam 
   """
}

process VAR_CAL {
  publishDir "${params.outdir}/variants", mode:'copy', overwrite:true

   input:
   tuple val(id), path(bam), path(ignore_index)
   path chr19

   output:
   tuple val(id), path("${id}.vcf")

   script:
   """
   module load bcftools

   bcftools mpileup -Ou -f ${chr19}/chr19.fa $bam | \
   bcftools call -mv -Ob -o ${id}.vcf
   """
}

