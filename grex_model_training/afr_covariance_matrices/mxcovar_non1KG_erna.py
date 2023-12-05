#!/usr/bin/env python

"""Python 3
"""

import argparse
import sys
import os.path
import sqlite3
import glob
import pandasvcf
import collections
import pandas
import numpy
import json

bad_locus_log = None

#Open the reference varID-rsID key file and convert to a Python dictionary
#dbsnp_df = pandas.read_csv("/home/bettimj/aldrich_rotation/mxcovar_keys/jti_varids_rsids_hg19_b150.txt", header = None, delim_whitespace = True)
#dbsnp_df.columns = ["varid", "rsid"]
#dbsnp_dict = dbsnp_df.set_index("rsid").T.to_dict()

def parse_varid(varid):
    words = varid.split("_")
    if len(words) == 1:
        words = varid.split(":")
    if len(words) < 4:
        sys.stderr.write(f"There was a problem with the format of varID - {varid}\n")
        sys.exit(1)
    (chr, pos, a1, a2) = words[0:4]
    pos = int(pos)
    return (chr, pos, a1, a2)

class Locus(object):
    def __init__(self, rsid, pos, ref, eff, weight):
        self.rsid = rsid
        self.pos = pos
        self.eff = eff
        self.ref = ref
        self.weight = weight
        self.genotypes = []

    def add_geno(self, ref, alt, geno):
        """We'll just ignore it if our alleles aren't aligned...for now"""
        alleles = [ref] + alt.split(",")

        try:
            al2 = str(alleles.index(self.eff))
            al1 = str(alleles.index(self.ref))
            lkup = {f"{al1}/{al1}":0, f"{al1}/{al2}":1, f"{al2}/{al1}":1, f"{al2}/{al2}":2}
        except:
            bad_locus_log.write("\t".join(["Invalid_Alleles", "%s/%s" % (self.ref, self.eff), self.rsid, "%s/%s" % (ref, alt)]))
            sys.stderr.write(f"Invalid Allele Combinations: {self.rsid}: {self.ref},{self.eff} -- {ref},{alt}\n")
            return
        cur_geno = []
        for gt in geno[5:]:
            if gt in lkup:
                cur_geno.append(lkup[gt])
            else:
                bad_locus_log.write("\t".join(["Oddball Genotype: ", f"{self.ref}/{self.eff}", self.rsid, f"{ref}/{alt}", gt]))
                # Force a hard return to skip for now
                return
        self.genotypes.append(cur_geno)

    def summarize_contents(self):
        print("\t" + "\t".join([
            self.rsid,
            str(self.pos),
            self.ref,
            self.eff,
            str(self.weight),
            str(len(self.genotypes))]))

# Directly ripped from Alvaro's code             -- https://github.com/hakyimlab/summary-gwas-imputation/blob/master/src/genomic_tools_lib/miscellaneous/matrices.py
def write_matrix_data(name, labels, matrix, file):
    """data is expected to be a list of (name, id_labels, matrix) tuples"""

    if len(labels) == 1:
        id = labels[0]

        file.write(f"{name} {id} {id} {str(matrix)}\n")
    else:
        for i in range(0, len(labels)):
            for j in range(i, len(labels)):
                value = matrix[i][j]
                id1 = labels[i]
                id2 = labels[j]
                file.write(f"{name} {id1} {id2} {str(value)}\n")

class Gene(object):
	def __init__(self, gene_name):
		self.name = gene_name
		self.df = None
		self.chrom = None

	def extract_weights(self, db):
		self.data = {}          # pos => Locus
		for (rsid, ref, eff, weight) in db.execute(f"SELECT rsid, ref_allele, eff_allele, weight FROM weights WHERE gene='{self.name}'"):
			print(rsid, ref, eff, weight)
#Append all rsIDs from the database to an array
			#varid = dbsnp_dict.get(rsid)
			#varid = varid["varid"]
			#print(varid)
			(chr, pos, a1, a2) = parse_varid(rsid)

            # If we have sex chromosomes or mt, we'll drop case
			chr = chr.replace("chr", "").lower()

			if self.chrom is None:
				self.chrom = chr
			if int(pos) in self.data:
				sys.stderr.write("Well, there is something already there: ", self.name, chr, pos)
				sys.exit(1)
			self.data[int(pos)] = Locus(rsid, int(pos), ref, eff, float(weight))

	def store_geno(self, pos, ref, alt, geno):
		if pos in self.data:
			self.data[pos].add_geno(ref, alt, geno)

	def summarize_contents(self):
		if len(self.data) > 1:
			print(f"{self.name} - {self.chrom}")
			for pos in sorted(self.data.keys()):
				self.data[pos].summarize_contents()

	def write_covariance(self, file):
		id_labels = []
		genotypes = []

		results = None
		if len(self.data) > 0:
			for pos in self.data:
				loc = self.data[pos]
				if len(loc.genotypes) > 0:
					id_labels.append(loc.rsid)
					genotypes += loc.genotypes
			cov_matrix = numpy.cov(numpy.array(genotypes))
			print(self.name, id_labels, numpy.array(genotypes), cov_matrix)
			write_matrix_data(self.name, id_labels, cov_matrix, file)


def LoadDbWeights(filename):
    db = sqlite3.connect(filename)
    cur = db.cursor()

    genes = {}
    for (genename,) in cur.execute("SELECT DISTINCT gene FROM weights"):
        genes[genename] = Gene(genename)

    # Chromosomes will be string
    chromosomes = collections.defaultdict(dict)
    for g in genes.keys():
        the_gene = genes[g]
        the_gene.extract_weights(cur)
        chrom = the_gene.chrom
        for pos in the_gene.data.keys():
            if pos not in chromosomes[chrom]:
                chromosomes[chrom][pos] = []
            chromosomes[chrom][pos].append(the_gene)

    return genes, chromosomes


def ParseArgs(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser(description='Process VCF file(s) and build covariance matrices ready for use with MetaXcan')
    parser.add_argument("--vcf", type=str, required=True, help='Single VCF file')
    parser.add_argument("--pop", '-p', type=argparse.FileType('r'), help='Text file with IDs to be kept for the final matrices')
    parser.add_argument("db", type=argparse.FileType('r'), nargs='+', help='One or more MetaXcan DB files to be used locus filtering')

    args = parser.parse_args()

    return args

def GetPopulation(args):
    population_ids = []
    for line in args.pop:
        # Catch IID, FID, etc
        if line[1:3].lower() != 'id':
            population_ids.append(line.strip().split()[0])
    return population_ids

def ListVCFs(args):
    "for now, let's just assume it's a list of VCFs"
    vcf_list = []
    if os.path.isdir(args.vcf):
        vcf_list = glob.glob(f"{args.vcf}/*.vcf") + glob.glob(f"{args.vcf}/*.vcf.gz")
    else:
        if os.path.exists(args.vcf):
            ext = os.path.splitext(args.vcf)[-1].lower()
            print(ext)
            if ext in ['.txt']:
                vcf_list = [x.strip() for x in args.vcf]
            elif ext in ['.vcf','.gz']:
                vcf_list = [args.vcf]
            else:
                sys.stderr.write(f"Unrecognized VCF file format, {args.vcf}")
        else:
            sys.stderr.write(f"Warning, {args.vcf} doesn't appear to be a directory nor a file!")
            sys.exit(1)
    return vcf_list

args = ParseArgs()
pop = GetPopulation(args)
vcfs = ListVCFs(args)

for db in args.db:
    prefix = db.name.split("/")[-1].split(".")[0]
    print(prefix)
    print(db.name)
    bad_locus_log = open("%s-bad_locus.log" % (prefix), 'w')
    covariance_file = open("%s.cov" % (prefix), 'w')
    covariance_file.write("GENE RSID1 RSID2 VALUE\n")
    genes, chromosomes = LoadDbWeights(db.name)

    for vcf in vcfs:
        vcf_data = pandasvcf.VCF(vcf, sample_id=pop)

        while vcf_data.get_vcf_df_chunk() == 0:
            for index,col in vcf_data.df.iterrows():
                chrom = col['CHROM']
                pos = col['POS']
                ref = col['REF']
                alt = col['ALT']
                if pos in chromosomes[chrom]:
                    for gene in chromosomes[chrom][pos]:
                        gene.store_geno(pos, ref, alt, col)
    visited = set() # Genes that have been printed, in case we have something weird
    for chr in sorted(chromosomes.keys()):
        for pos in sorted(chromosomes[chr]):
            for gene in chromosomes[chr][pos]:
                if gene.name not in visited:
                    gene.write_covariance(covariance_file)
                    visited.add(gene.name)
