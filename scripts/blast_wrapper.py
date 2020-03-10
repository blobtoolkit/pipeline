#!/usr/bin/env python3
"""Run BLAST with long sequences split into chunks."""

import math
import sys
import time
import shlex
import logging
import re
import glob
import os.path
from collections import defaultdict
from random import shuffle
from itertools import groupby
from multiprocessing import Pool
from subprocess import Popen, PIPE, run, CalledProcessError

docs = """
Run BLAST

Usage: ./blast_wrapper.py -query FASTA -db BLASTDB

Options:
    -program blastn           BLAST program to use.
    -db BLASTDB               BLAST database.
    -query FASTAFILE          query sequence FASTA file.
    -taxidlist TAXIDFILE      list of taxids to include (requires v5 BLAST DB).
    -multiprocessing          use python multiprocessing for num_threads.
    -chunk 100000             sequences greater than CHUNK bp will be split.
    -overlap 500              length of overlap when splitting sequences.
    -max-chunks 10            maximum number of chunks to split a sequence into.
    -chunks FASTAFILE.chunks  chunked sequence filename.
    -num-threads 16           number of threads.
    -evalue 1e-25             BLAST evalue.
    -max-target-seqs 10       BLAST max_target_seqs.
    -raw FASTAFILE.out.raw    raw output filename.
    -out FASTAFILE.out        BLAST output filename (processed from raw output).
    -nohit FASTAFILE.nohit    query sequences with no hit to BLAST DB.
    -window_masker_db FILE    windowmasker optimized binary format counts file.
    ...                       any additional blast parameters.
"""

logger_config = {
    'level': logging.INFO,
    'format': '%(asctime)s [%(levelname)s] line %(lineno)d %(message)s',
    'filemode': 'w'
}
try:
    logger_config.update({'filename': snakemake.log[0]})
except NameError as err:
    pass
logging.basicConfig(**logger_config)
logger = logging.getLogger()


def parse_args():
    """Parse command line arguments."""
    # using custom parser to allow passthrough of unspecified arguments
    script_params = {
        '-program': 'blastn',
        '-query': None,
        '-chunk': 100000,
        '-multiprocessing': 'False',
        '-overlap': 500,
        '-max_chunks': 10,
        '-num_threads': 16,
        '-raw': None,
        '-out': '.out',
        '-nohit': '.nohit',
        '-chunks': None
    }
    blast_params = {
        '-db': None,
        '-evalue': '1e-25',
        '-max_target_seqs': '10',
        # '-outfmt': '6 qseqid staxids bitscore std',
        '-max_hsps': '1'
    }
    try:
        script_params['-chunks'] = snakemake.output.chunks
        script_params['-query'] = snakemake.input.fasta
        blast_params['-db'] = "%s/%s" % (snakemake.params.dir, snakemake.wildcards.name)
        blast_params['-taxidlist'] = snakemake.input.taxids
        script_params['-multiprocessing'] = str(snakemake.params.multiprocessing)
        script_params['-chunk'] = int(snakemake.params.chunk)
        script_params['-overlap'] = int(snakemake.params.overlap)
        script_params['-max_chunks'] = int(snakemake.params.max_chunks)
        script_params['-num_threads'] = int(snakemake.threads)
        blast_params['-evalue'] = str(snakemake.params.evalue)
        blast_params['-max_target_seqs'] = str(snakemake.params.max_target_seqs)
        try:
            script_params['-raw'] = snakemake.output.raw
        except AttributeError:
            pass
        try:
            script_params['-nohit'] = snakemake.output.nohit
        except AttributeError:
            pass
        try:
            script_params['-out'] = snakemake.output.out
        except AttributeError:
            pass
        script_params['-max_target_seqs'] = blast_params['-max_target_seqs']
        if script_params['-multiprocessing'] == 'False':
            blast_params['-num_threads'] = str(script_params['-num_threads'])
        blast_list = [item for k in blast_params for item in (k, blast_params[k])]
        return (script_params, blast_list)
    except AttributeError:
        script_params['-query'] = snakemake.input.fasta
        blast_params['-db'] = "%s/%s" % (snakemake.params.dir, snakemake.wildcards.name)
        blast_params['-taxidlist'] = snakemake.input.taxids
        script_params['-multiprocessing'] = str(snakemake.params.multiprocessing)
        script_params['-num_threads'] = int(snakemake.threads)
        blast_params['-evalue'] = str(snakemake.params.evalue)
        blast_params['-max_target_seqs'] = str(snakemake.params.max_target_seqs)
        try:
            script_params['-raw'] = snakemake.output.raw
        except AttributeError:
            pass
        try:
            script_params['-nohit'] = snakemake.output.nohit
        except AttributeError:
            pass
        try:
            script_params['-out'] = snakemake.output.out
        except AttributeError:
            pass
        script_params['-max_target_seqs'] = blast_params['-max_target_seqs']
        if script_params['-multiprocessing'] == 'False':
            blast_params['-num_threads'] = str(script_params['-num_threads'])
        blast_list = [item for k in blast_params for item in (k, blast_params[k])]
        return (script_params, blast_list)
    except NameError as err:
        logger.info(err)
        logger.info('Parsing parameters from command line')
    try:
        pair = False
        for arg in sys.argv:
            if re.match(r'-\w', arg):
                pair = arg
            elif pair:
                if pair in script_params:
                    script_params[pair] = arg
                else:
                    blast_params[pair] = arg
                pair = False
        if script_params['-multiprocessing'] == 'False':
            blast_params['-num_threads'] = str(script_params['-num_threads'])
        for arg in ['-raw', '-out', '-nohit']:
            if script_params[arg] and re.match('.', script_params[arg]):
                script_params[arg] = script_params['-query']+script_params[arg]
        blast_list = [item for k in blast_params for item in (k, blast_params[k])]
        script_params['-max_target_seqs'] = blast_params['-max_target_seqs']
        for arg in ['-chunk', '-overlap', '-max_chunks', '-num_threads', '-max_target_seqs']:
            script_params[arg] = int(script_params[arg])
    except Exception as err:
        logger.error(err)
        logger.error('Unable to parse parameters')
        print(docs)
        exit(1)
    if not os.path.exists(script_params['-query']):
        logger.error("'%s' is not a valid query file" % script_params['-query'])
        print(docs)
        exit(1)
    if '-window_masker_db' in blast_params and not os.path.exists(blast_params['-window_masker_db']):
        logger.error("'%s' is not a valid windowmasker file" % blast_params['-window_masker_db'])
        print(docs)
        exit(1)
    if not glob.glob("%s.*" % blast_params['-db']):
        logger.error("'%s' is not a valid database file" % blast_params['-db'])
        print(docs)
        exit(1)
    return (script_params, blast_list)


def chunk_size(value):
    """Calculate nice value for chunk size."""
    mag = math.floor(math.log10(value))
    first = int(str(value)[:2]) + 1
    chunk_size = first * pow(10, mag-1)
    return chunk_size


def split_list(input, size):
    """Yield successive subsets from list."""
    for i in range(0, len(input), size):
        yield input[i:i + size]


def chunk_fasta(fastafile, chunk=math.inf, overlap=0, max_chunks=math.inf):
    """Read FASTA file one sequence at a time and split long sequences into chunks."""
    cmd = "cat %s" % fastafile
    # TODO: read gzipped files if needed
    # cmd = "pigz -dc %s" % fastafile
    title = ''
    seq = ''
    segment = chunk + overlap
    with Popen(shlex.split(cmd), encoding='utf-8', stdout=PIPE, bufsize=4096) as proc:
        faiter = (x[1] for x in groupby(proc.stdout, lambda line: line[0] == '>'))
        for header in faiter:
            title = header.__next__()[1:].strip().split()[0]
            seq = ''.join(map(lambda s: s.strip(), faiter.__next__()))
            seq_length = len(seq)
            my_chunk = chunk
            if seq_length > segment:
                n = (seq_length + chunk) // chunk
                if n > max_chunks:
                    my_chunk = chunk_size(seq_length / max_chunks)
                    n = max_chunks
                for i in range(0, seq_length, my_chunk):
                    subseq = seq[i:i+segment]
                    yield {'title': title, 'seq': subseq, 'chunks': n, 'start': i}
            else:
                yield {'title': title, 'seq': seq, 'chunks': 1, 'start': 0}


def read_fasta(fastafile):
    """Read FASTA file one sequence at a time."""
    cmd = "cat %s" % fastafile
    # TODO: read gzipped files if needed
    # cmd = "pigz -dc %s" % fastafile
    title = ''
    seq = ''
    with Popen(shlex.split(cmd), encoding='utf-8', stdout=PIPE, bufsize=4096) as proc:
        faiter = (x[1] for x in groupby(proc.stdout, lambda line: line[0] == '>'))
        for header in faiter:
            title = header.__next__()[1:].strip().split()[0]
            seq = ''.join(map(lambda s: s.strip(), faiter.__next__()))
            seq_length = len(seq)
            yield {'title': title, 'seq': seq, 'chunks': 1, 'start': 0}


def run_blast(seqs, cmd, blast_list, index, batches):
    """Run blast on seqs."""
    logger.info("running BLAST on %d sequences in batch %d of %d" % (len(seqs), index, batches))
    try:
        input = ''
        for seq in seqs:
            input += ">%s_-_%d\n" % (seq['title'], seq['start'])
            input += "%s\n" % seq['seq']
        cmd += ' -lcase_masking -outfmt "6 qseqid staxids bitscore std"'
        logger.info(shlex.split(cmd)+blast_list)
        p = run(shlex.split(cmd)+blast_list, stdout=PIPE, stderr=PIPE, input=input, encoding='ascii')
        return p
    except Exception as err:
        logger.error("Unable to run %s" % cmd)
        logger.info(blast_list)
        logger.error(err)
        exit(1)


def parse_raw_output(output, outfile, count):
    """Process raw output to standard BLAST output format."""
    try:
        lines = defaultdict(dict)
        chunk_counts = defaultdict(int)
        for line in output:
            fields = line.split('\t')
            if fields[0]:
                name, start = re.split('_-_', fields[0])
                fields[0] = name
                fields[3] = name
                fields[9] = str(int(fields[9])+int(start))
                fields[10] = str(int(fields[10])+int(start))
                if start not in lines[name]:
                    lines[name][start] = []
                    chunk_counts[name] += 1
                lines[name][start].append('\t'.join(fields))
    except Exception as err:
        logger.error(err)
        logger.error("Unable to parse raw output")
        exit(1)
    try:
        with open(outfile, 'w') as ofh:
            for name in chunk_counts.keys():
                length = len(lines[name])
                n = (int(count)+length-1) // length
                for start in lines[name].keys():
                    for i in range(n):
                        if i < len(lines[name][start]):
                            ofh.write("%s\n" % str(lines[name][start][i]))
    except Exception as err:
        logger.error(err)
        logger.error("Unable to write output to %s" % outfile)
        exit(1)


if __name__ == '__main__':
    (script_params, blast_list) = parse_args()
    logger.info("Running blast_wrapper")
    try:
        seqs = []
        names = set()
        if script_params['-chunk']:
            for seq in chunk_fasta(script_params['-query'],
                                   chunk=script_params['-chunk'],
                                   overlap=script_params['-overlap'],
                                   max_chunks=script_params['-max_chunks']):
                if not re.match('^N+$', seq['seq']):
                    names.add(seq['title'])
                    seqs.append((seq))
            if script_params['-chunks']:
                chunked = ''
                for seq in seqs:
                    chunked += ">%s_-_%d\n" % (seq['title'], seq['start'])
                    chunked += "%s\n" % seq['seq']
                with open(script_params['-chunks'], 'w') as ofh:
                    ofh.writelines(chunked)
        else:
            for seq in read_fasta(script_params['-query']):
                if not re.match('^N+$', seq['seq']):
                    names.add(seq['title'])
                    seqs.append((seq))
        n_chunks = len(seqs)
        if n_chunks == 0:
            logger.error("no sequences found in query file '%s'" % script_params['-query'])
            exit(1)
        shuffle(seqs)
        subset_length = math.ceil(n_chunks / script_params['-num_threads'])
        min_length = subset_length / script_params['-num_threads']
        while subset_length > script_params['-num_threads'] and subset_length > min_length:
            subset_length //= 2
        if script_params['-multiprocessing'] != 'False':
            pool = Pool(script_params['-num_threads'])
            jobs = []
            pool_error = None
        output = []

        def close_pool(error):
            """Close pool on error."""
            if error:
                logger.error(error)
                pool.terminate()

        def blast_callback(p):
            """Process BLAST chunk."""
            global output
            logger.info('entering callback')
            try:
                p.check_returncode()
            except CalledProcessError as err:
                close_pool(err)
                logger.error(p.stderr)
                logger.error("Unable to run %s" % script_params['-program'])
                exit(1)
            result = ''
            if p.stderr:
                logger.info(p.stderr)
            lines = []
            if p.stdout:
                lines = p.stdout.strip('\n').split('\n')
            logger.info("Finished processing batch with %d BLAST hits" % len(lines))
            for line in lines:
                fields = line.split('\t')
                if fields[0]:
                    output.append('\t'.join(fields))

        index = 1
        batches = math.ceil(len(seqs) / subset_length)
        try:
            if script_params['-multiprocessing'] == 'False':
                logger.info('running single BLAST')
                p = run_blast(seqs, script_params['-program'], blast_list, 1, 1)
                logger.info('calling callback')
                blast_callback(p)
            else:
                for subset in split_list(seqs, subset_length):
                    proc = pool.apply_async(run_blast, (subset, script_params['-program'], blast_list, index, batches), callback=blast_callback)
                    jobs.append(proc)
                    index += 1
                pool.close()
                pool.join()
                for job in jobs:
                    job.get()
                    job.wait()
        except Exception as err:
            logger.error(err)
            logger.error("Unable to process chunks")
            exit(1)
        logger.info("Finished BLASTing %d chunks" % (index - 1))
        if script_params['-raw']:
            logger.info("Writing raw output to file '%s'" % script_params['-raw'])
            try:
                with open(script_params['-raw'], 'w') as ofh:
                    ofh.writelines('\n'.join(output))
                for line in output:
                    name = line.split('_-_')[0]
                    if name in names:
                        names.remove(name)
            except Exception as err:
                logger.error(err)
                logger.error("Unable to write raw output to %s" % script_params['-raw'])
                exit(1)
        if script_params['-out']:
            logger.info("Writing output to file '%s'" % script_params['-out'])
            parse_raw_output(output, script_params['-out'], script_params['-max_target_seqs'])
        if script_params['-nohit']:
            logger.info("Writing nohit IDs to file '%s'" % script_params['-nohit'])
            try:
                with open(script_params['-nohit'], 'w') as ofh:
                    ofh.writelines('\n'.join(names))
            except Exception as err:
                logger.error(err)
                logger.error("Unable to write nohit IDs to %s" % script_params['-nohit'])
                exit(1)
    except Exception as err:
        logger.error(err)
        exit(1)
