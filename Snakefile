import re
from pprint import pprint
from glob import glob

configfile: "config.json"
localrules: all, mkdir

sample_id_2_gz_pair = {}

def get_sample_ids(input_dir):
    global sample_id_2_gz_pair
    for x in glob('{}/*'.format(input_dir)):
        if not x.startswith('.'):
            m = re.match(r'input/(.+?)_S\d+_L001_R[12]_001.fastq.gz', x)
            if not m:
                raise Exception("unexpected file name: {}".format(x))
            sample_id = m.group(1)
            if sample_id not in sample_id_2_gz_pair:
                sample_id_2_gz_pair[sample_id] = [x]
            else:
                sample_id_2_gz_pair[sample_id].append(x)
    for k, v in sample_id_2_gz_pair.items():
        if len(v) != 2:
            raise Exception('{} not paired'.format(k))
        v.sort()
    return sorted(list(sample_id_2_gz_pair.keys()))

sample_ids = get_sample_ids('input')

rule all:
    input:
        "table/raw_OTU_table_uncollapsed.txt",
        "table/raw_OTU_table_collapsed.txt",
        "table/QC_table.txt"

rule mkdir:
    output: touch("mkdir.done")
    shell: "mkdir -p unzip join filter fq_2_fa strip collapse blast blast_parse table"

rule unzip:
    input:
        lambda wildcards: sample_id_2_gz_pair[wildcards.sample_id]
    output:
        "unzip/{sample_id}_R1.fq",
        "unzip/{sample_id}_R2.fq",
    shell:
        "mkdir -p unzip && gunzip -c {input[0]} > {output[0]} && gunzip -c {input[1]} > {output[1]}"

rule join:
    input: rules.unzip.output
    output: "join/{sample_id}.assembled.fastq"
    version: config["tool_versions"]["pear"]
    shell: "tools/pear/{version}/pear -f {input[0]} -r {input[1]} -o join/{wildcards.sample_id}"

rule filter:
    input: rules.join.output
    output: "filter/{sample_id}.fq"
    version: config["tool_versions"]["fastx"]
    params:
        percentage = config["parameters"]["quality_filtering"]["percentage"],
        qscore = config["parameters"]["quality_filtering"]["qscore"]
    shell: "if [ ! -s {input} ]; then touch {output}; else tools/fastx/{version}/fastq_quality_filter -i {input} -o {output} -q {params.qscore} -p {params.percentage} -Q33 -v; fi"

rule fq_2_fa:
    input: rules.filter.output
    output: "fq_2_fa/{sample_id}.fa"
    version: config["tool_versions"]["fastx"]
    shell: "if [ ! -s {input} ]; then touch {output}; else tools/fastx/{version}/fastq_to_fasta -i {input} -o {output} -n -r -v -Q33; fi"

rule strip:
    input: rules.fq_2_fa.output
    output: "strip/{sample_id}.fa"
    params:
        forward_primer = config["parameters"]["strip"]["forward_primer"],
        reverse_primer = config["parameters"]["strip"]["reverse_primer"]
    shell: "if [ ! -s {input} ]; then touch {output}; else python scripts/strip.py -f {params.forward_primer} -r {params.reverse_primer} -i {input} -o {output}; fi"

rule collapse:
    input: rules.strip.output
    output: "collapse/{sample_id}.fa"
    version: config["tool_versions"]["fastx"]
    shell: "if [ ! -s {input} ]; then touch {output}; else tools/fastx/{version}/fastx_collapser -i {input} -o {output} -v; fi"

rule mkdb:
    input: expand("database/{name}/{version}/{file}", name=config["database"]["name"], version=config["database"]["version"], file=config["database"]["fa_file"])
    output: touch("mkdb.done")
    version: config["tool_versions"]["blast+"]
    params: config["database"]["name"]
    shell: "tools/blast+/{version}/bin/makeblastdb -in {input} -out {params} -dbtype 'nucl' -input_type fasta"

rule blast:
    input:
        query = rules.collapse.output,
        mkblastdb_done = "mkdb.done"
    output: "blast/{sample_id}.txt"
    version: config["tool_versions"]["blast+"]
    params:
        identity = config["parameters"]["blast"]["identity"],
        coverage = config["parameters"]["blast"]["coverage"],
        db_name = config["database"]["name"]
    shell: 'tools/blast+/{version}/bin/blastn -query {input.query} -task megablast -db {params.db_name} -max_target_seqs 5000 -perc_identity {params.identity} -qcov_hsp_perc {params.coverage} -evalue 1e-6 -outfmt "7 qacc sacc qstart qend sstart send length pident qcovhsp qcovs" -out {output}'

rule blast_parse:
    input: "blast{sample_id}.txt"
    output: "blast_parse{sample_id}.txt"
    shell: "python scripts/blast_parse.py -i {input} -o {output}"

rule make_otu_table:
    input:
        blast_res = expand("blast_parse/{sample_id}.txt", sample_id=sample_ids),
        taxonomy = expand("database/{name}/{version}/{file}", name=config["database"]["name"], version=config["database"]["version"], file=config["database"]["tax_file"])
    output:
        uncollapsed_table = "table/raw_OTU_table_uncollapsed.txt",
        collapsed_table = "table/raw_OTU_table_collapsed.txt"
    shell: "python scripts/make_OTU_table.py -b {input.blast_res} -t {input.taxonomy} -u {output.uncollapsed_table} -c {output.collapsed_table}"

rule qc:
    input:
        raw_count = expand("unzip/{sample_id}_R1.fq", sample_id=sample_ids),
        join_count = expand("join/{sample_id}.assembled.fastq", sample_id=sample_ids),
        filter_count = expand("filter/{sample_id}.fq", sample_id=sample_ids),
        blast_count = expand("blast_parse/{sample_id}.txt", sample_id=sample_ids)
    output: "table/QC_table.txt"
    shell: "python scripts/qc.py -r {input.raw_count} -j {input.join_count} -f {input.filter_count} -b {input.blast_count} -o {output}"

rule clean:
    params: config["database"]["name"]
    shell: "rm -rf {params}.* mkdir.done mkdb.done unzip join filter fq_2_fa strip collapse blast blast_parse table slurm-*.out"

