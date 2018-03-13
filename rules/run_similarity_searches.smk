rule run_blastn:
    """
    Run NCBI blastn to search nucleotide database with assembly query.
    """
    input:
        fna='{assembly}.fna',
        db='{name}.root.{root}{masked}.nal'
    output:
        '{assembly}.blastn.{name}.root.{root}{masked}.out'
    wildcard_constraints:
        root='\d+',
        masked='[minus\d\.]+'
    params:
        db=lambda wc: "%s.root.%s%s" % (wc.name,wc.root,wc.masked),
        evalue=lambda wc:similarity[wc.name]['evalue'],
        max_target_seqs=lambda wc:similarity[wc.name]['max_target_seqs']
    threads: 32
    priority: 10
    shell:
        '{ENV} \
        blastn \
            -query {input.fna} \
            -db {params.db} \
            -outfmt "6 qseqid staxids bitscore std" \
            -max_target_seqs {params.max_target_seqs} \
            -max_hsps 1 \
            -evalue {params.evalue} \
            -num_threads {threads} \
            > {output}'

rule run_blastx:
    """
    Run NCBI blastx to search protein database with assembly query.
    """
    input:
        fna='{assembly}.fna',
        db='{name}.root.{root}{masked}.pal'
    output:
        '{assembly}.blastx.{name}.root.{root}{masked}.out'
    wildcard_constraints:
        root='\d+',
        masked='[minus\d\.]+'
    params:
        db=lambda wc: "%s.root.%s%s" % (wc.name,wc.root,wc.masked),
        evalue=lambda wc:similarity[wc.name]['evalue'],
        max_target_seqs=lambda wc:similarity[wc.name]['max_target_seqs']
    threads: 32
    priority: 20
    shell:
        '{ENV} \
        blastx \
            -query {input.fna} \
            -db {params.db} \
            -outfmt "6 qseqid staxids bitscore std" \
            -max_target_seqs {params.max_target_seqs} \
            -max_hsps 1 \
            -evalue {params.evalue} \
            -num_threads {threads} \
            > {output}'


rule run_diamond_blastx:
    """
    Run Diamond blastx to search protein database with assembly query.
    """
    input:
        fna='{assembly}.fna',
        db='{name}.root.{root}{masked}.dmnd'
    output:
        '{assembly}.diamond.{name}.root.{root}{masked}.out'
    wildcard_constraints:
        root='\d+',
        masked='[minus\d\.]+'
    params:
        db=lambda wc: "%s.root.%s%s" % (wc.name,wc.root,wc.masked),
        evalue=lambda wc:similarity[wc.name]['evalue'],
        max_target_seqs=lambda wc:similarity[wc.name]['max_target_seqs']
    threads: 32
    shell:
        '{ENV} \
        diamond blastx \
            --query {input.fna} \
            --db {params.db} \
            --outfmt 6 \
            --sensitive \
            --max-target-seqs {params.max_target_seqs} \
            --evalue {params.evalue} \
            --threads {threads} \
            > {output}'

rule run_blobtools_taxify:
    """
    Add taxonomy information to Diamond similarity search results.
    """
    input:
        dmnd='{assembly}.diamond.{name}.root.{root}{masked}.out',
        idmap=lambda wc: "%s/full/%s.taxid_map.gz" % (similarity[wc.name]['local'],wc.name),
        nodes="%s/data/nodesDB.txt" % config['settings']['blobtools']
    output:
        '{assembly}.diamond.{name}.root.{root}{masked}.taxified.out'
    params:
        idmap=lambda wc: "%s/full/%s.taxid_map" % (similarity[wc.name]['local'],wc.name),
    wildcard_constraints:
        root='\d+',
        masked='[minus\d\.]+'
    conda:
        '../envs/blobtools.yaml'
    threads: 1
    shell:
        '{ENV} \
        pigz -dc {input.idmap} > {params.idmap} && \
        PATH=' + config['settings']['blobtools'] + ':$PATH && \
        blobtools taxify \
            -f {input.dmnd} \
            -m {params.idmap} \
            -s 0 \
            -t 1'
