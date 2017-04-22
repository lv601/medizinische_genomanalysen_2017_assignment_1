#! /usr/bin/env python3

import mysql.connector
import os
import pybedtools
import pysam
import subprocess
import pybedtools
import pickle
import pprint

__author__ = 'mgollobich'


class Assignment1:
    def __init__(self):
        self.target = "HMBS"
        self.gene = self.fetch_gene_coordinates(self,"hg19",)
        self.geneinfo = ['BDNF', 'NM_001258209', '11', 118955586, 118964259, '+', 13,
                         "b'118955586,118958964,118959344,118959791,118959926,118960392,118960699,118960899,118962122,118962834,118963467,118963644,118963819"]
        self.bamfile = os.path.join(os.getcwd(), "HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam")
        self.bedtoolsfile = pybedtools.BedTool('HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam')
        self.coveragefile = "Coveragefile"

        for line in open('MYFILE.TXT'):
            if line.startswith("(u'HMBS'"):
                self.info = line
                print (line)


        # Check if bam file exist, if not download
        if not os.path.isfile(self.bamfile):
            subprocess.call(["wget", "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/"
                                     "alignment/HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120"
                                     "522.bam"])
        # Check if bai file exist, if not create
        if not os.path.isfile(self.bamfile + ".bai"):
            subprocess.call(["samtools", "index", self.bamfile])

        self.samfile = pysam.AlignmentFile(self.bamfile, "rb")
        self.reads = list(self.samfile.fetch(self.gene.chrom, self.gene.txStart, self.gene.txEnd))

    @staticmethod
    def fetch_gene_coordinates(self, genome_reference):

        print("Connecting to UCSC to fetch data of target gene\t", self.target)
        cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu', user='genomep', passwd='password',
                                      db=genome_reference)
        cursor = cnx.cursor()

        query_fields = ["refGene.name2",
                        "refGene.name",
                        "refGene.chrom",
                        "refGene.txStart",
                        "refGene.txEnd",
                        "refGene.strand",
                        "refGene.exonCount",
                        "refGene.exonStarts",
                        "refGene.exonEnds"]

        query = "SELECT DISTINCT %s from refGene" % ",".join(query_fields) + \
                " WHERE refGene.name2=" + '"' + self.target + '"' ""
        cursor.execute(query)

        for row in cursor:
            self.gene = gene(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8])

        cursor.close()
        cnx.close()

        print("Done fetching data\n")

        return self.gene

    def get_sam_header(self):
        print("get sam header:")
        for key, value in self.samfile.header['HD'].items():
            if key == "SO":
                print("Sorting order of alignments (SO): ", value)
            if key == "VN":
                print("Format version (VN): ", value)
            if key == "GO":
                print("Grouping of alignments (GO): ", value)

    def get_properly_paired_reads_of_gene(self):
        print('The proper paired reads are:')
        n = 0
        for read in self.samfile.fetch("11", self.geneinfo[3], self.geneinfo[4]):
            if read.is_proper_pair:
                print("-properly_paired_reads_of_gene:")
                n += 1
                print(read)
        print("Number of properly paired reads: {}".format(n))

    def get_gene_reads_with_indels(self):
        print("Genes with Indels: \n")
        liste=[]
        for read in self.samfile:
            columns = str(read).split("\t")
            if "I" in str(columns[5]) or "D" in str(columns[5]):
                #print("-gene_reads_with_indels:")
                liste.append(columns[5])
        # Find indels in unmapped reads by using Cigar
        i = 0
        for read in self.reads:
            if not read.is_unmapped:
                cigar = read.cigar
                for (type, length) in cigar:
                    # Insertion or Deletion
                    if (type == 1) or (type == 2):
                        i += 1
        if i == 0:
            print("No gene reads with indesls found \n")
        else:
            print("Number of gene reads with indels:")
            print(i," ... \n")
            print(liste[:10])

    def calculate_total_average_coverage(self):
        print("Start calculating total average coverage...\n")
        coverage = self.bedtoolsfile.genome_coverage(bg=True)
        #if not os.path.isfile(self.coveragefile):
        #    coverage = self.bedtoolsfile.genome_coverage(bg=True)
        #    pickle.dump(coverage, open(self.coveragefile, "wb"))
        #else:
        #    coverage = pickle.load(open(self.coveragefile, "rb"))

        #print(coverage[:10])
        i = 0
        #for i, line in enumerate(coverage):
            #print(coverage)

        i = 0
        average = 0

        for line in coverage:
            # Converting from str to float
            number = float(line[3])
            average += number
            i += 1

        coverage = average/i

        print("Total average coverage:")
        print(coverage, "... \n")

    def calculate_gene_average_coverage(self):


        print("Start calculating gene average coverage...\n")
        coverage = self.bedtoolsfile.genome_coverage(bg=True)

        average = 0
        i = 0

        for line in coverage:
            number = float(line[3])
            cbeg = int(line[1])

            if cbeg > self.gene.txStart:
                if int(line[2]) <= self.gene.txEnd:
                    average += number
                    i += 1

        coverage = average / i

        print("Total gene average coverage:")
        print(coverage, '... \n')

    def get_number_mapped_reads(self):
        i = 0
        for read in self.reads:
            if not read.is_unmapped:
                i += 1
        if i == 0:
            print("no mapped reads reads found")
        else:
            print("Number of mapped reads: ",i,"\n")


    def get_gene_symbol(self):
        print("Gene symbol: ",self.target)

    def get_region_of_gene(self):
        print('Region of gene - Start: {} End: {} Lenght: {}'.format(self.geneinfo[3],self.geneinfo[4],self.geneinfo[4]-self.geneinfo[3]))
        #print('Start:',self.geneinfo[3],'End: ',self.geneinfo[4],'Lenght: ',self.geneinfo[4]-self.geneinfo[3])

    def get_number_of_exons(self):
        print("Number of exons: ",self.geneinfo[6])


    def print_summary(self):
        print("Summary - Gene of Interest: ", self.target, "\n")
        self.get_sam_header()
        self.get_properly_paired_reads_of_gene()
        self.get_gene_reads_with_indels()
        self.calculate_total_average_coverage()
        self.calculate_gene_average_coverage()
        self.get_number_mapped_reads()
        self.get_gene_symbol()
        self.get_region_of_gene()
        self.get_number_of_exons()
        print("\n",'\n')

class gene:
    def __init__(self, name2, name, chrom, txStart, txEnd, strand, exonCount, exonStarts, exonEnds):
        self.name2 = name2
        self.name = name
        self.chrom = chrom[10:]
        self.txStart = txStart
        self.txEnd = txEnd
        self.strand = strand
        self.exonCount = exonCount
        self.exonStarts = str(exonStarts).lstrip("b'").rstrip(",'").split(",")
        self.exonEnds = str(exonEnds).lstrip("b'").rstrip(",'").split(",")



if __name__ == '__main__':
    print("Assignment 1: ", __author__)
    assignment1 = Assignment1()
    assignment1.print_summary()
