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
    tuple val(chrX), val(sex_id), path(vcf), path(vcf_index)

    output:
    tuple stdout, val(chrX), val(sex_id), path(vcf), path(vcf_index)

    """
    tabix -l ${vcf} 
    """
}

process rm_chr_name_ref{
    cache "lenient"
    cpus 1
    memory "4GB"
    time "00:30:00"

    scratch true
    input:
    tuple val(chr_name), val(chrX), val(sex_id), path(vcf), path(vcf_index)

    output:
    tuple val(chr_name), val(chrX), val(sex_id), path("*.vcf.gz"), path("*.vcf.gz.tbi")
    
    script:
    if (chrX == false) {
    """
    echo "${chr_name}" > chroms1.txt
    chr=${chr_name}
    echo "\${chr:3}" > chroms2.txt

    paste chroms1.txt chroms2.txt > chr_name_conv.txt  

    bcftools annotate --rename-chrs chr_name_conv.txt $vcf -Oz -o ${vcf.getBaseName()}.vcf.gz
    bcftools index --tbi ${vcf.getBaseName()}.vcf.gz
    """
    } else {
    """
    echo "${chr_name}" > chroms1.txt
    echo "X" > chroms2.txt

    paste chroms1.txt chroms2.txt > chr_name_conv.txt  

    bcftools annotate --rename-chrs chr_name_conv.txt $vcf -Oz -o ${vcf.getBaseName()}.vcf.gz
    bcftools index --tbi ${vcf.getBaseName()}.vcf.gz
    """
    }
}


process get_study_chr_names{
    cache "lenient"
    cpus 1
    memory "4GB"
    time "00:30:00"
    scratch true

    input:
    tuple val(chrX), val(sex_id), path(vcf), path(vcf_index)

    output:
    tuple stdout, val(chrX), val(sex_id), path(vcf), path(vcf_index)

    """
    tabix -l ${vcf}
    """
}

process rm_chr_name_study{
    cache "lenient"
    cpus 1
    memory "4GB"
    time "00:30:00"
    scratch true  

    input:
    tuple val(chr_name), val(chrX), val(sex_id), path(vcf), path(vcf_index)

    output:
    tuple val(chr_name), val(chrX), val(sex_id), path("*.vcf.gz"), path("*.vcf.gz.tbi")
    
    script:
    if (chrX == false) {
    """
    echo "${chr_name}" > chroms1.txt
    chr=${chr_name}
    echo "\${chr:3}" > chroms2.txt

    paste chroms1.txt chroms2.txt > chr_name_conv.txt  

    bcftools annotate --rename-chrs chr_name_conv.txt $vcf -Oz -o ${vcf.getBaseName()}.vcf.gz
    bcftools index --tbi ${vcf.getBaseName()}.vcf.gz
    """
    } else {
    """
    echo "${chr_name}" > chroms1.txt
    echo "X" > chroms2.txt

    paste chroms1.txt chroms2.txt > chr_name_conv.txt  

    bcftools annotate --rename-chrs chr_name_conv.txt $vcf -Oz -o ${vcf.getBaseName()}.vcf.gz
    bcftools index --tbi ${vcf.getBaseName()}.vcf.gz
    """
    }
}

process get_chunks_study {
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

process get_chunks_ref {
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
process convert_ref_vcf{
    //errorStrategy 'retry'
    //maxRetries 3
    cache "lenient"
    cpus 4
    memory "32GB"
    time "5h"
    scratch true

    input:
    tuple val(chr_name), val(chrX), val(sex_id), path(chunk), path(ref_vcf), path(ref_vcf_index)
    
    output:
    tuple val(chr_name), val(chrX), val(sex_id), path(chunk), path( "*.m3vcf.gz")
    
    publishDir "minimac_m3vcfs/", pattern: "*.m3vcf.gz", mode: "copy"

    """
    ${params.minimac3} --refHaps $ref_vcf --cpus 4 --processReference --prefix ${ref_vcf.getBaseName()}
    """
}

process minimac_imputation{
    //errorStrategy 'retry'
    //maxRetries 3
    cache "lenient"
    cpus 4
    memory "32GB"
    time "5h"
    scratch true  

    input:
    tuple val(chr_name), val(chrX), val(sex_id), file(ref_vcf), file(study_vcf), file(study_vcf_index)
     
    output:
    tuple path("*.vcf.gz*"), path("*.info")

    publishDir "imputed_vcfs/", pattern: "*.vcf.gz*", mode: "copy"
    publishDir "imputed_info/", pattern: "*.info", mode: "copy"

    script:
    if (chrX == false) {
    """
    ${params.minimac4} --refHaps $ref_vcf --haps $study_vcf  --prefix ${study_vcf.getBaseName()} --meta --ignoreDuplicates
    echo "${chr_name}" > chroms1.txt
    chr=${chr_name}
    echo "\${chr:3}" > chroms2.txt

    paste chroms2.txt chroms1.txt > chr_name_conv.txt  

    bcftools annotate --rename-chrs chr_name_conv.txt ${study_vcf.getBaseName()}.dose.vcf.gz -Oz -o ${study_vcf.getBaseName()}.imp.dose.vcf.gz
    bcftools index --tbi ${study_vcf.getBaseName()}.imp.dose.vcf.gz

    bcftools annotate --rename-chrs chr_name_conv.txt ${study_vcf.getBaseName()}.empiricalDose.vcf.gz -Oz -o ${study_vcf.getBaseName()}.imp.empiricalDose.vcf.gz
    bcftools index --tbi ${study_vcf.getBaseName()}.imp.empiricalDose.vcf.gz
    """
    } else {
    """
    ${params.minimac4} --refHaps $ref_vcf --haps $study_vcf --cpus 4 --prefix ${study_vcf.getBaseName()} --meta --ignoreDuplicates

    echo "chrX" > chroms1.txt
    echo "X" > chroms2.txt

    paste chroms2.txt chroms1.txt > chr_name_conv.txt  

    bcftools annotate --rename-chrs chr_name_conv.txt ${study_vcf.getBaseName()}.dose.vcf.gz -Oz -o ${study_vcf.getBaseName()}.imp.dose.vcf.gz
    bcftools index --tbi ${study_vcf.getBaseName()}.imp.dose.vcf.gz

    bcftools annotate --rename-chrs chr_name_conv.txt ${study_vcf.getBaseName()}.empiricalDose.vcf.gz -Oz -o ${study_vcf.getBaseName()}.imp.empiricalDose.vcf.gz
    bcftools index --tbi ${study_vcf.getBaseName()}.imp.empiricalDose.vcf.gz
    """
    }
}
workflow {

        ref_ch = Channel.fromPath(params.ref_vcf_path).map{ vcf -> [ vcf.name.toString().tokenize('.').contains('chrX'), vcf.name.toString().tokenize('.').contains('female'), vcf, vcf + ".tbi" ] }
        study_ch = Channel.fromPath(params.study_vcf_path).map{ vcf -> [ vcf.name.toString().tokenize('.').contains('chrX_array'), vcf.name.toString().tokenize('.').contains('female'), vcf, vcf + ".tbi" ] }


        ref_vcfs = get_ref_chr_names(ref_ch)
        study_vcfs = get_study_chr_names(study_ch)
        
        ref_rm_chr_vcfs = rm_chr_name_ref(ref_vcfs)
        study_rm_chr_vcfs = rm_chr_name_study(study_vcfs)

        ref_cnv_vcfs = convert_ref_vcf(ref_rm_chr_vcfs)

        imputation_ch = ref_cnv_vcfs.join(study_rm_chr_vcfs, by:[0, 1, 2])

        minimac_imputation(imputation_ch)
}