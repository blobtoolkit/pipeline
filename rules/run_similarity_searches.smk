rule chunk_fasta:
    """
    split fasta file into chunks.
    """
    input:
        '{assembly}.fasta'
    output:
        '{assembly}.fasta.chunked'
    params:
        chunk=config['settings']['blast_chunk'],
        overlap=config['settings']['blast_overlap'],
    conda:
         '../envs/py3.yaml'
    threads: 1
    resources:
        threads=1
    script:
        '../scripts/chunk_fasta.py'

rule run_blastn:
    """
    Run NCBI blastn to search nucleotide database with assembly query.
    """
    input:
        fasta='{assembly}.fasta',
        db=lambda wc: "%s/%s.nal" % (similarity[wc.name]['local'],wc.name),
        taxids='{name}.root.{root}{masked}.taxids'
    output:
        raw='{assembly}.blastn.{name}.root.{root}{masked}.out.raw',
        nohit=temp('{assembly}.blastn.{name}.root.{root}{masked}.nohit')
    wildcard_constraints:
        root='\d+',
        masked='.[fm][ulins\d\.]+'
    params:
        db=lambda wc: "%s/%s" % (similarity[wc.name]['local'],wc.name),
        evalue=lambda wc:similarity[wc.name]['evalue'],
        max_target_seqs=lambda wc:similarity[wc.name]['max_target_seqs'],
        path=config['settings']['blast_path'],
        chunk=config['settings']['blast_chunk'],
        overlap=config['settings']['blast_overlap'],
        max_chunks=config['settings']['blast_max_chunks']
    conda:
         '../envs/pyblast.yaml'
    threads: 64
    resources:
        threads=64
    script:
        '../scripts/blast_wrapper.py'
    # shell:
    #     '{params.path}/blastn \
    #         -query {input.fasta} \
    #         -db {params.db} \
    #         -outfmt "6 qseqid staxids bitscore std" \
    #         -max_target_seqs {params.max_target_seqs} \
    #         -max_hsps 1 \
    #         -evalue {params.evalue} \
    #         -num_threads {threads} \
    #         -taxidlist {input.taxids} \
    #         > {output}'

rule unchunk_blast_results:
    """
    reformat blast results from chunked input.
    """
    input:
        '{assembly}.{algorithm}.{name}.root.{root}{masked}.out.raw'
    output:
        '{assembly}.{algorithm}.{name}.root.{root}{masked}.out'
    params:
        max_target_seqs=lambda wc:similarity[wc.name]['max_target_seqs']
    conda:
         '../envs/py3.yaml'
    threads: 1
    resources:
        threads=1
    script:
        '../scripts/unchunk_blast.py'

# rule run_blastn:
#     """
#     Run NCBI blastn to search nucleotide database with assembly query.
#     """
#     input:
#         fasta='{assembly}.fasta',
#         db=lambda wc: "%s/%s.nal" % (similarity[wc.name]['local'],wc.name),
#         taxids='{name}.root.{root}{masked}.taxids'
#     output:
#         '{assembly}.blastn.{name}.root.{root}{masked}.out'
#     wildcard_constraints:
#         root='\d+',
#         masked='.[fm][ulins\d\.]+'
#     params:
#         db=lambda wc: "%s/%s" % (similarity[wc.name]['local'],wc.name),
#         evalue=lambda wc:similarity[wc.name]['evalue'],
#         max_target_seqs=lambda wc:similarity[wc.name]['max_target_seqs'],
#         blast_chunk=config['settings']['blast_chunk'],
#         blast_overlap=config['settings']['blast_overlap'],
#         path=config['settings']['blast_path']
#     conda:
#          '../envs/pyblast.yaml'
#     threads: 32
#     resources:
#         threads=32
#     script:
#         '../scripts/blast_wrapper.py'

rule extract_nohit_sequences:
    """
    Run seqtk to extract nohit sequences from assembly.
    """
    input:
        fasta='{assembly}.fasta',
        nohit='{assembly}.blastn.{name}.root.{root}{masked}.nohit'
    output:
        temp('{assembly}.blastn.{name}.root.{root}{masked}.fasta.nohit')
    conda:
         '../envs/blast.yaml'
    threads: 1
    resources:
        threads=1
    shell:
        'seqtk subseq {input.fasta} {input.nohit} > {output}'

rule run_blastx:
    """
    Run NCBI blastx to search protein database with assembly query.
    """
    input:
        fasta='{assembly}.blastn.nt.root.{root}{masked}.fasta.nohit',
        db='blast/{name}.root.{root}{masked}.pal'
    output:
        '{assembly}.blastx.{name}.root.{root}{masked}.out'
    wildcard_constraints:
        root='\d+',
        masked='.[fm][ulins\d\.]+'
    params:
        db=lambda wc: "%s.root.%s%s" % (wc.name,wc.root,wc.masked),
        evalue=lambda wc:similarity[wc.name]['evalue'],
        max_target_seqs=lambda wc:similarity[wc.name]['max_target_seqs'],
        path=config['settings']['blast_path']
    conda:
         '../envs/blast.yaml'
    threads: 64
    resources:
        threads=64
    shell:
        'cd blast && \
        {params.path}/blastx \
            -query ../{input.fasta} \
            -db {params.db} \
            -outfmt "6 qseqid staxids bitscore std" \
            -max_target_seqs {params.max_target_seqs} \
            -max_hsps 1 \
            -evalue {params.evalue} \
            -num_threads {threads} \
            > ../{output}'

rule run_diamond_blastx:
    """
    Run Diamond blastx to search protein database with assembly query.
    """
    input:
        fasta='{assembly}.blastn.nt.root.{root}{masked}.fasta.nohit',
        db='{name}.root.{root}{masked}.dmnd'
    output:
        '{assembly}.diamond.{name}.root.{root}{masked}.out'
    wildcard_constraints:
        root='\d+',
        masked='.[fm][ulins\d\.]+',
        assembly='\w+'
    params:
        db=lambda wc: "%s.root.%s%s" % (wc.name,wc.root,wc.masked),
        evalue=lambda wc:similarity[wc.name]['evalue'],
        max_target_seqs=lambda wc:similarity[wc.name]['max_target_seqs']
    conda:
         '../envs/diamond.yaml'
    threads: 64
    resources:
        threads=64
    shell:
        'diamond blastx \
            --query {input.fasta} \
            --db {params.db} \
            --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
            --sensitive \
            --max-target-seqs {params.max_target_seqs} \
            --evalue {params.evalue} \
            --threads {threads} \
            > {output}'

# rule run_blobtools_taxify:
#     """
#     Add taxonomy information to Diamond similarity search results.
#     """
#     input:
#         dmnd='{assembly}.diamond.{name}.root.{root}{masked}.out',
#         split=lambda wc: "%s/split/%s.done" % (similarity[wc.name]['local'],wc.name),
#         lists='blast/{name}.root.{root}{masked}.lists'
#     output:
#         '{assembly}.diamond.{name}.root.{root}{masked}.taxified.out'
#     params:
#         indir=lambda wc: "%s/split/%s" % (similarity[wc.name]['local'],wc.name),
#         idmap=lambda wc: "%s/%s.taxid_map" % (config['settings']['tmp'],wc.name),
#         path=config['settings']['blobtools_path']
#     wildcard_constraints:
#         root='\d+',
#         masked='.[fm][ulins\d\.]+'
#     conda:
#         '../envs/blobtools.yaml'
#     threads: 1
#     resources:
#         threads=1
#     shell:
#         'parallel --no-notice -j {threads} \
#             "gunzip -c {params.indir}/{{}}.taxid_map.gz" \
#             :::: {input.lists} > {params.idmap} && \
#         {params.path}/blobtools taxify \
#             -f {input.dmnd} \
#             -m {params.idmap} \
#             -s 0 \
#             -t 1 && \
#         rm {params.idmap}'
