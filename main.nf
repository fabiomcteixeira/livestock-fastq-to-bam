nextflow.enable.dsl=2
import groovy.json.JsonSlurper

process NORMALIZE_LB {
    tag "$animal/$id"
    cpus 1
    memory '4 GB'
    input:
    tuple val(animal), val(id), path(r1), path(r2)
    output:
    tuple val(animal), val(id), path("${id}.normalized_R1.fastq.gz"), path("${id}.normalized_R2.fastq.gz")
    script:
    """
    normalize_lb.py '${r1}' '${r2}' '${id}.normalized_R1.fastq.gz' '${id}.normalized_R2.fastq.gz' > normalization.tsv
    """
}

process MAP {
    tag "$animal/$id"
    cpus 16
    memory '48 GB'
    publishDir { "${params.outdir}/04_mapping/${id}" }, mode: 'symlink', overwrite: true, failOnError: true
    input:
    tuple val(animal), val(id), path(r1), path(r2), val(library)
    path reference_files
    output:
    tuple val(animal), val(id), path("${id}.sorted.bam"), path("${id}.sorted.bam.bai")
    path 'bwa.log', emit: log
    path 'versions.txt', emit: versions
    script:
    def ref = file(params.reference).name
    """
    ("\$PILOT_BWA" 2>&1 || true) > versions.txt
    "\$PILOT_SAMTOOLS" --version >> versions.txt
    "\$PILOT_BWA" mem -t 10 -R '@RG\\tID:${id}\\tSM:${animal}\\tLB:${library}\\tPL:ILLUMINA' \
      '${ref}' '${r1}' '${r2}' 2> bwa.log |
      "\$PILOT_SAMTOOLS" sort -@ 4 -m 2G -o '${id}.sorted.bam' -
    "\$PILOT_SAMTOOLS" quickcheck -v '${id}.sorted.bam'
    "\$PILOT_SAMTOOLS" index -@ 4 '${id}.sorted.bam'
    """
}

process MERGE_RUNS {
    tag "$id"
    cpus 4
    memory '16 GB'
    publishDir { "${params.outdir}/04_merged/${id}" }, mode: 'symlink', overwrite: true, failOnError: true
    input:
    tuple val(id), path(bams)
    output:
    tuple val(id), path("${id}.merged.bam"), path("${id}.merged.bam.bai")
    script:
    def ordered = bams instanceof List ? bams.sort { it.name } : [bams]
    def inputs = ordered.collect { "'${it}'" }.join(' ')
    def merge = ordered.size() == 1 ? "ln -s ${inputs} '${id}.merged.bam'" : "samtools merge -@ 3 -o '${id}.merged.bam' ${inputs}"
    """
    ${merge}
    samtools quickcheck -v '${id}.merged.bam'
    samtools index -@ 3 '${id}.merged.bam'
    """
}

process DEDUP {
    tag "$id"
    cpus 4
    memory '64 GB'
    publishDir { "${params.outdir}/05_duplicates/${id}" }, mode: 'symlink', overwrite: true, failOnError: true
    input:
    tuple val(id), path(bam), path(bai)
    output:
    tuple val(id), path("${id}.dedup.bam"), path("${id}.dedup.bam.bai"), emit: bam
    path '*.metrics.txt', emit: metrics
    path 'picard.log', emit: log
    path 'versions.txt', emit: versions
    script:
    """
    java -jar '${params.picard_jar}' MarkDuplicates --version > versions.txt 2>&1 || {
        grep -Eq '^Version:[[:space:]]*[0-9]' versions.txt || { cat versions.txt >&2; exit 1; }
    }
    mkdir -p tmp
    java -Xmx24g -XX:ActiveProcessorCount=4 -jar '${params.picard_jar}' MarkDuplicates I='${bam}' \
      O='${id}.dedup.bam' M='${id}.duplicates.metrics.txt' \
      REMOVE_DUPLICATES=false REMOVE_SEQUENCING_DUPLICATES=true CREATE_INDEX=false \
      VALIDATION_STRINGENCY=LENIENT TMP_DIR=tmp 2> picard.log
    samtools quickcheck -v '${id}.dedup.bam'
    samtools index -@ 3 '${id}.dedup.bam'
    """
}

process TRIM {
    tag "$animal/$id"
    cpus 8
    memory '16 GB'
    publishDir { "${params.outdir}/02_trim/${id}" }, mode: 'symlink', overwrite: true, failOnError: true
    input:
    tuple val(animal), val(id), path(r1), path(r2), path(adapters), val(quality_options)
    output:
    tuple val(animal), val(id), path("${id}_R1P.fq.gz"), path("${id}_R2P.fq.gz"), emit: paired
    path '*U.fq.gz', emit: unpaired
    tuple val(id), path("${id}.trim.json"), emit: counts
    path 'trimming.log', emit: log
    path 'versions.txt', emit: versions
    script:
    """
    java -jar '${params.trimmomatic_jar}' -version > versions.txt 2>&1
    java -Xmx12g -jar '${params.trimmomatic_jar}' PE -threads ${task.cpus} -phred33 '${r1}' '${r2}' \
      '${id}_R1P.fq.gz' '${id}_R1U.fq.gz' \
      '${id}_R2P.fq.gz' '${id}_R2U.fq.gz' \
      ILLUMINACLIP:${adapters}:2:30:10:1:true ${quality_options} 2> trimming.log
    qc_tools.py trim '${animal}' '${id}' trimming.log > ${id}.trim.json
    """
}

process FASTQC_RAW {
    tag "$id"
    cpus 2
    memory '4 GB'
    publishDir { "${params.outdir}/01_fastqc_raw/${id}" }, mode: 'symlink', overwrite: true
    input:
    tuple val(animal), val(id), path(r1), path(r2)
    output:
    tuple val(animal), val(id), path(r1), path(r2), path("${id}.raw_qc.json"), emit: checked
    path '*_fastqc.*', emit: reports
    path 'versions.txt', emit: versions
    script:
    """
    '${params.fastqc}' --version > versions.txt
    '${params.fastqc}' --threads 2 --outdir . '${r1}' '${r2}'
    qc_tools.py qc '${animal}' '${id}' raw *_fastqc.zip > '${id}.raw_qc.json'
    """
}

process FASTQC_TRIM {
    tag "$id"
    cpus 2
    memory '4 GB'
    publishDir { "${params.outdir}/03_fastqc_trimmed/${id}" }, mode: 'symlink', overwrite: true
    input:
    tuple val(animal), val(id), path(r1), path(r2)
    output:
    tuple val(animal), val(id), path(r1), path(r2), path("${id}.trim_qc.json"), emit: checked
    path '*_fastqc.*', emit: reports
    path 'versions.txt', emit: versions
    script:
    """
    '${params.fastqc}' --version > versions.txt
    '${params.fastqc}' --threads 2 --outdir . '${r1}' '${r2}'
    qc_tools.py qc '${animal}' '${id}' trim *_fastqc.zip > '${id}.trim_qc.json'
    """
}

process BAM_STATS {
    tag "$id/$kind"
    cpus 2
    memory '4 GB'
    publishDir { "${params.outdir}/stats/${kind}/${id}" }, mode: 'symlink', overwrite: true
    input:
    tuple val(animal), val(id), val(kind), path(bam), path(bai)
    output:
    path "${id}.${kind}.json", emit: summary
    path '*.txt', emit: details
    path '*.tsv', emit: depth
    script:
    """
    samtools flagstat '${bam}' > flagstat.txt
    samtools stats '${bam}' > samtools_stats.txt
    samtools view -c -F 2308 '${bam}' > mapped.txt
    samtools depth -aa -s -q 20 -Q 20 -G 3844 '${bam}' | depth_summary.py > depth.tsv
    qc_tools.py bam '${animal}' '${id}' '${kind}' mapped.txt depth.tsv > '${id}.${kind}.json'
    """
}

process REPORT {
    executor 'local'
    cpus 1
    memory '1 GB'
    publishDir "${params.outdir}/tables", mode: 'copy', overwrite: true
    input:
    path records
    path manifest
    output:
    path '*.tsv'
    script:
    """
    qc_tools.py report '${manifest}' *.json
    """
}

workflow {
    if (!(params.stage in ['raw_qc','trim_qc','full'])) error 'Unknown stage'
    rows = file(params.manifest, checkIfExists:true).splitCsv(header:true, sep:'\t')
    expected = rows.countBy { it.animal_id }
    approved = file(params.approvals, checkIfExists:true).readLines().collect { it.trim() }.findAll { it && !it.startsWith('#') } as Set
    unknown = approved - (rows.collect { it.run } as Set)
    if (unknown) error "Unknown approved runs: ${unknown}"
    overrides = file(params.trim_overrides, checkIfExists:true).splitCsv(header:true, sep:'\t').collectEntries { [(it.run):it.options] }
    if (overrides.keySet().any { !(it in rows.collect { r -> r.run }) }) error 'Unknown run in trim_overrides.tsv'
    raw = Channel.fromList(rows).map { r -> tuple(r.animal_id, r.run, file(r.R1), file(r.R2)) }
    libraries = rows.collectEntries { [(it.run): (it.library_id ?: it.run)] }
    formats = rows.collectEntries { [(it.run):it.header_format] }
    FASTQC_RAW(raw)
    reports = FASTQC_RAW.out.checked.map { a,id,r1,r2,j -> j }
    if (params.stage != 'raw_qc') {
        eligible = FASTQC_RAW.out.checked.filter { a,id,r1,r2,j ->
            def q = new JsonSlurper().parseText(j.text)
            q.adapter_pass || (approved.contains(id) && q.valid)
        }.map { a,id,r1,r2,j -> tuple(a,id,r1,r2) }
        NORMALIZE_LB(eligible.filter { a,id,r1,r2 -> formats[id]=='hash_mate' })
        ready = eligible.filter { a,id,r1,r2 -> formats[id]=='standard' }.mix(NORMALIZE_LB.out)
        // A per-run value avoids re-trimming unaffected runs when approvals change.
        trim_inputs = ready.map { a,id,r1,r2 ->
            def quality = overrides[id] ?: (params.quality_options)
            def adapters = approved.contains(id) ? file(params.combined_adapters) : file(params.adapters)
            tuple(a,id,r1,r2,adapters,quality)
        }
        TRIM(trim_inputs)
        FASTQC_TRIM(TRIM.out.paired)
        reports = reports.mix(TRIM.out.counts.map { id,j -> j }, FASTQC_TRIM.out.checked.map { a,id,r1,r2,j -> j })
        if (params.stage == 'full') {
            passed = FASTQC_TRIM.out.checked.filter { a,id,r1,r2,j ->
                def q = new JsonSlurper().parseText(j.text)
                q.adapter_pass && q.quality_pass && q.valid && q.total_reads > 0
            }.map { a,id,r1,r2,j -> tuple(a,id,r1,r2,libraries[id]) }
            refs = ['', '.amb','.ann','.bwt','.pac','.sa'].collect { file(params.reference+it, checkIfExists:true) }
            MAP(passed, refs)
            // Emit each animal as soon as ALL expected run BAMs exist.
            grouped = MAP.out[0].map { a,id,bam,bai -> tuple(groupKey(a,expected[a]),bam) }
                .groupTuple().map { key,bams -> tuple(key.getGroupTarget(),bams) }
            MERGE_RUNS(grouped)
            DEDUP(MERGE_RUNS.out)
            stats_inputs = MAP.out[0].map { a,id,bam,bai -> tuple(a,id,'run',bam,bai) }
                .mix(MERGE_RUNS.out.map { a,bam,bai -> tuple(a,a,'merged',bam,bai) },
                     DEDUP.out.bam.map { a,bam,bai -> tuple(a,a,'final',bam,bai) })
            BAM_STATS(stats_inputs)
            reports = reports.mix(BAM_STATS.out.summary)
        }
    }
    REPORT(reports.collect(), file(params.manifest))
}
