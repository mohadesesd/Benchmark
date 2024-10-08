#!/usr/bin/env nextflow
/*
* AUTHOR: Mohadese Sayahian Dehkordi, <mohadese.sayahiandehkordi@mail.mcgill.ca>
* Adapted and modified get_chunk and ligate from https://github.com/CERC-Genomic-Medicine/shapeit4_pipeline/blob/main/Phasing.nf by Daniel Taliun, PhD <daniel.taliun@mcgill.ca>
* VERSION: 3.0
* YEAR: 2023
*/


process get_ref_chr_names{
    cache "lenient"
    cpus 1
    memory "4GB"
    time "00:30:00"
    scratch true
    input:
    tuple val(sex_id), path(vcf), path(vcf_index)

    output:
    tuple stdout, val(sex_id), path(vcf), path(vcf_index)

    """
    chrom=`bcftools index -s ${vcf} | cut -f1`
	printf "\${chrom}"
    """
}


process get_study_chr_names{
    cache "lenient"
    cpus 1
    memory "4GB"
    time "00:30:00"
    scratch true
    input:
    tuple val(sex_id), path(vcf), path(vcf_index)

    output:
    tuple stdout, val(sex_id), path(vcf), path(vcf_index)

    """
    chrom=`bcftools index -s ${vcf} | cut -f1`
	printf "\${chrom}"
    """
}

process get_chunks {
	//executor "local"
	cache "lenient"
	cpus 1
	memory "4GB"
    time "00:30:00"

	input:
	tuple val(chromosome), val(sex_id), path(vcf), path(vcf_index)

    output:
    tuple path("*.chunk"), val(chromosome), val(sex_id), path(vcf), path(vcf_index)

	"""
    chrom=`bcftools index -s ${vcf} | cut -f1`
    start_bp=`bcftools view -HG ${vcf} | head -n1 | cut -f2`
	stop_bp=`bcftools index -s ${vcf} | cut -f2`
       
    extend=0
	for i in `seq \${start_bp} ${params.window} \${stop_bp}`; do
		if [ \${extend} -eq 0 ]; then
			chunk_start=\$((${params.flank} > i ? 0 : i - ${params.flank}))
		fi
		chunk_stop=\$((i + ${params.window} + ${params.flank}))
		n=`bcftools view -HG ${vcf} \${chrom}:\${chunk_start}-\${chunk_stop} | wc -l`
		if [ \${n} -gt 0 ]; then
			printf "\${chrom}\t\${chunk_start}\t\${chunk_stop}\t\${n}\n" > \${chrom}_\${chunk_start}_\${chunk_stop}.chunk
			extend=0
		else
			extend=1
		fi
	done
	if [ \${extend} -eq 1 ]; then
		printf "\${chrom}\t\${chunk_start}\t\${chunk_stop}\t\${n}\n" > \${chrom}_\${chunk_start}_\${chunk_stop}.chunk
	fi
	"""
}


process eagle_phasing{
    //errorStrategy 'retry'
    //maxRetries 3
    cache "lenient"
    cpus 4
    memory "16GB"
    time "05:00:00"
    scratch true
    input:
    tuple val(chromosome), val(sex_id), path(ref_vcf), path(ref_vcf_index), path(chunk), path(study_vcf), path(study_vcf_index) 
        
    output:
    tuple val(chromosome), val(sex_id), path("*.phased.final.vcf.gz"), path("*.phased.final.vcf.gz.tbi")
    publishDir "phased_vcfs_chunks/", pattern: "*.vcf.gz*", mode: "copy"

    script:
    if(params.chromosomeX == true){
    """
    echo "${chromosome}" > chroms1.txt
    echo "chrX" > chroms2.txt

    paste chroms1.txt chroms2.txt > chr_name_conv.txt   
    
    bcftools annotate --rename-chrs chr_name_conv.txt $study_vcf -Oz -o ${study_vcf.getBaseName()}.vcf.gz
    bcftools index --tbi ${study_vcf.getBaseName()}.vcf.gz

    bcftools annotate --rename-chrs chr_name_conv.txt $ref_vcf -Oz -o ${ref_vcf.getBaseName()}.vcf.gz
    bcftools index --tbi ${ref_vcf.getBaseName()}.vcf.gz


    chrom=`head -n1 ${chunk} | cut -f1`
    start_bp=`head -n1 ${chunk} | cut -f2`
	stop_bp=`head -n1 ${chunk} | cut -f3`
    stop_w_overlap=\$((stop_bp+5000000))
    bcftools view -r \${chrom}:\${start_bp}-\${stop_bp} ${study_vcf.getBaseName()}.vcf.gz -Oz -o study.\${chrom}_\${start_bp}_\${stop_bp}.vcf.gz 
    bcftools index --tbi study.\${chrom}_\${start_bp}_\${stop_bp}.vcf.gz
    bcftools view -r \${chrom}:\${start_bp}-\${stop_w_overlap} ${ref_vcf.getBaseName()}.vcf.gz -Oz -o ref.\${chrom}_\${start_bp}_\${stop_w_overlap}.vcf.gz 
    bcftools index --tbi ref.\${chrom}_\${start_bp}_\${stop_w_overlap}.vcf.gz
    ${params.eagle} --vcfRef ref.\${chrom}_\${start_bp}_\${stop_w_overlap}.vcf.gz  --vcfTarget study.\${chrom}_\${start_bp}_\${stop_bp}.vcf.gz --numThreads 4 --geneticMapFile ${params.genetic_map} --outPrefix \${chrom}_\${start_bp}_\${stop_bp}.phased --allowRefAltSwap --vcfOutFormat z
    bcftools index --tbi \${chrom}_\${start_bp}_\${stop_bp}.phased.vcf.gz
    paste chroms2.txt chroms1.txt > chr_name_conv.txt   

    bcftools annotate --rename-chrs chr_name_conv.txt \${chrom}_\${start_bp}_\${stop_bp}.phased.vcf.gz -Oz -o \${chrom}_\${start_bp}_\${stop_bp}.phased.final.vcf.gz
    bcftools index --tbi \${chrom}_\${start_bp}_\${stop_bp}.phased.final.vcf.gz

    """
    } else {
    """
    chrom=`head -n1 ${chunk} | cut -f1`
    start_bp=`head -n1 ${chunk} | cut -f2`
	stop_bp=`head -n1 ${chunk} | cut -f3`
    stop_w_overlap=\$((stop_bp + 5000000))
    bcftools view -r \${chrom}:\${start_bp}-\${stop_bp} ${study_vcf} -Oz -o study.\${chrom}_\${start_bp}_\${stop_bp}.vcf.gz 
    bcftools index --tbi study.\${chrom}_\${start_bp}_\${stop_bp}.vcf.gz
    bcftools view -r \${chrom}:\${start_bp}-\${stop_w_overlap} ${ref_vcf} -Oz -o ref.\${chrom}_\${start_bp}_\${stop_w_overlap}.vcf.gz 
    bcftools index --tbi ref.\${chrom}_\${start_bp}_\${stop_w_overlap}.vcf.gz
    ${params.eagle} --vcfRef ref.\${chrom}_\${start_bp}_\${stop_w_overlap}.vcf.gz  --vcfTarget study.\${chrom}_\${start_bp}_\${stop_bp}.vcf.gz --numThreads 4 --geneticMapFile ${params.genetic_map} --outPrefix \${chrom}_\${start_bp}_\${stop_bp}.phased.final --allowRefAltSwap --vcfOutFormat z
    bcftools index --tbi \${chrom}_\${start_bp}_\${stop_bp}.phased.final.vcf.gz
    """
    }

}
process ligate {
	
	cache "lenient"
	scratch true
	cpus 1
	memory "32G"
	time "12h"

	input:
	tuple val(chromosome), val(sex_id), path(phased_vcfs), path(phased_vcfs_index)

	output:
	tuple path("*.phased.vcf.gz"), path("*.phased.vcf.gz.tbi")

	publishDir "phased_vcfs/", pattern: "*.vcf.gz*", mode: "copy"

    script:
    if(params.chromosomeX == true){
	"""
	for f in ${phased_vcfs}; do echo \${f}; done | sort -V > files_list.txt
	bcftools concat -f files_list.txt -l -Oz -o ${chromosome}.phased.vcf.gz
	bcftools index --tbi ${chromosome}.phased.vcf.gz
	"""
    }else {
    """
    if [ "${sex_id}" == "true" ]; then
        sex="female"
    else
        sex="male"
    fi
    for f in ${phased_vcfs}; do echo \${f}; done | sort -V > files_list.txt
	bcftools concat -f files_list.txt -l -Oz -o ${chromosome}.\${sex}.phased.vcf.gz
	bcftools index --tbi ${chromosome}.\${sex}.phased.vcf.gz
	"""
    }
}

workflow {
        ref_ch = Channel.fromPath(params.ref_vcf_path).map{ vcf -> [ vcf.name.toString().tokenize('.').contains('female'), vcf, vcf + ".tbi" ] }
        study_ch = Channel.fromPath(params.study_vcf_path).map{ vcf -> [ vcf.name.toString().tokenize('.').contains('female'), vcf, vcf + ".tbi" ] }

        ref_vcfs = get_ref_chr_names(ref_ch)
        study_vcfs = get_study_chr_names(study_ch)
        study_chunks = get_chunks(study_vcfs)
        chunks_all = study_chunks.flatMap { chunks, chromosome, sex_id, vcf, vcf_index ->
        chunks.collect { chunk -> [chromosome, sex_id, chunk, vcf, vcf_index] }
        }
        phasing_ch = ref_vcfs.combine(chunks_all, by:[0, 1])
        phased_chunks = eagle_phasing(phasing_ch)
        ligate(phased_chunks.groupTuple(by:[0, 1]))
}
