#! /usr/bin/env python3

import argparse
import os

import misc.utils as utils

default_sjdb_overhang = 50
default_num_procs = 2
default_star_executable = 'STAR'
default_mem = "31G"

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='This script creates the STAR reference necessary for the genomic mapping.')
    parser.add_argument('sjdb', help='The splice junction database file (gtf format)')
    parser.add_argument('fasta', help='Fasta file containing the reference genome')
    parser.add_argument('out', help='The path to the indexed genome')
    parser.add_argument('-p', '--num-procs', help="Number of threads to use", type=int, 
        default=default_num_procs)
    parser.add_argument('-m', '--mem', help="The amount of memory to use", default=default_mem)
    parser.add_argument('--overhang', help='The overhang to use for splice junctions '
        '(STAR --sjdbOverhang)', type=int, default=default_sjdb_overhang)
    parser.add_argument('--star-executable', help="The name of the STAR executable",
        default=default_star_executable)
    args = parser.parse_args()

    if not os.path.exists(args.out):
        os.makedirs(args.out)

    # check if the proper STAR executable is available
    utils.check_programs_exist([args.star_executable])

    mem = utils.human2bytes(args.mem)
    cmd = "{} --runMode genomeGenerate --genomeDir {} --genomeFastaFiles {} --sjdbGTFfile {} --sjdbOverhang {} --runThreadN {} --limitGenomeGenerateRAM {}".format(args.star_executable, args.out, args.fasta, args.sjdb, args.overhang, args.num_procs, mem)
    utils.check_call(cmd)

if __name__ == '__main__':
    main()




    #/software/STAR/STAR_2.4.0h/STAR --runMode genomeGenerate --genomeDir . --genomeFastaFiles c_elegans.PRJNA13758.WS246.genomic_softmasked.fa --sjdbGTFfile /data/Genomes/C_elegans/Caenorhabditis_elegans.WBcel235.    76.combined.gtf --sjdbOverhang 50 --runThreadN 8
