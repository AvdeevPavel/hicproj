import pysam

from gfa_graph.asm_graph import AsmGraph


def parse_cigar(cigar):
    matched_bp = 0
    algn_ref_span = 0
    algn_read_span = 0
    read_len = 0
    clip5_ref = 0
    clip3_ref = 0

    if cigar != '*':
        cur_num = 0
        for char in cigar:
            charval = ord(char)
            if charval >= 48 and charval <= 57:
                cur_num = cur_num * 10 + (charval - 48)
            else:
                if char == 'M':
                    matched_bp += cur_num
                    algn_ref_span += cur_num
                    algn_read_span += cur_num
                    read_len += cur_num
                elif char == 'I':
                    algn_read_span += cur_num
                    read_len += cur_num
                elif char == 'D':
                    algn_ref_span += cur_num
                elif char == 'S' or char == 'H':
                    read_len += cur_num
                    if matched_bp == 0:
                        clip5_ref = cur_num
                    else:
                        clip3_ref = cur_num

                cur_num = 0

    return {
        'clip5_ref': clip5_ref,
        'clip3_ref': clip3_ref,
        'cigar': cigar,
        'algn_ref_span': algn_ref_span,
        'algn_read_span': algn_read_span,
        'read_len': read_len,
        'matched_bp': matched_bp,
    }

UNMAPPED_POS = 0

# def parse_alternative_algns(samcols):
#     alt_algns = []
#     for col in samcols[11:]:
#         if not col.startswith('XA:Z:'):
#             continue
#
#         for SA in col[5:].split(';'):
#             if not SA:
#                 continue
#             SAcols = SA.split(',')
#
#             chrom = SAcols[0]
#             strand = '-' if SAcols[1]<0 else '+'
#
#             cigar = parse_cigar(SAcols[2])
#             NM = SAcols[3]
#
#             pos = UNMAPPED_POS
#             if strand == '+':
#                 pos = int(SAcols[1])
#             else:
#                 pos = int(SAcols[1]) + cigar['algn_ref_span']
#
#             alt_algns.append({
#                 'chrom': chrom,
#                 'pos': pos,
#                 'strand': strand,
#                 'mapq': mapq,
#                 'is_mapped': True,
#                 'is_unique': False,
#                 'is_linear': None,
#                 'cigar': cigar,
#                 'NM': NM,
#                 'dist_to_5': cigar['clip5_ref'] if strand == '+' else cigar['clip3_ref'],
#             })
#
#     return supp_algns
# 7938
# 628
# 8566
# 0
# 222
#
# 0
# 1256
# 0
# 0
# 7938 -- mapped
# 628 -- unmapped
# 1263
# 0 6950


'''
When a read matches in its entirety, with an equal score in multiple locations, one of the locations is picked at 
random, is labeled as primary, will be given a mapping quality of zero and will have an XA tag that contains the 
alternative locations (this is identical to how bwa aln worked)

When different, non-overlapping regions of a read align with high scores to different, non-linear locations in 
the genome, the higher score alignment will be labeled as primary, the others may be reported as secondary alignments. 
There is some threshold on how many of these secondary alignments will be reported 

When two complementary regions of a read (the two pieces add up to the full read) align to two different, non-linear 
genomic locations one of the alignment will be labeled as primary, the other as supplementary alignment 
'''


class BAMParser(object):

    def __init__(self, bam_file: str, graph: AsmGraph) -> None:
        self.path_to_bam = bam_file
        self.graph = graph

        # total number of alignment records
        self.total_number = 0
        # n = non-mapped, m = mapped, each letter is read from pair
        self.nn_count, self.nm_count, self.mm_count = 0, 0, 0
        # supp = supplementary, sec = secondary, prime = prim
        self.supp_count, self.sec_count, self.prime_count = 0, 0, 0
        # mapq0 is interesting since it may have perfect alignments but with XA tag, posq can be further filtered
        self.mapq0_count, self.posq_count = 0, 0

    def parse(self) -> None:
        bam_stream = pysam.AlignmentFile(self.path_to_bam, 'rb')

        for record in bam_stream:
            self.total_number += 1
            # if record.is_unmapped or record.mate_is_unmapped:
            if record.flag & 0x0004 and record.flag & 0x0008:
                # both reads are unmapped - we are not interested
                self.nn_count += 1
            elif (record.flag & 0x0004) != (record.flag & 0x0008):
                # one of two reads is unmapped - we are not interested
                self.nm_count += 1
            else:
                self.mm_count += 1
                if record.flag & 0x0800:  # record.is_supplementary:
                    # TODO small overlaps because of repeats or chimeric stuff
                    self.supp_count += 1
                elif record.flag & 0x0100:  # record.is_secondary:
                    # TODO may cause by overlaps in unitigs or complete overlaps
                    self.sec_count += 1
                else:
                    self.prime_count += 1
                    if record.mapping_quality == 0:
                        self.mapq0_count += 1
                    else:
                        self.posq_count += 1

    def write_statistics_info_into_file(self, out_file: str) -> None:
        with open(out_file, 'w') as out:
            out.write("Total number of alignment records = {0}\n".format(self.total_number))
            out.write("\tBoth reads in pair are unmapped {0}\n".format(self.nn_count))
            out.write("\tOne of reads in pair are unmapped {0}\n".format(self.nm_count))
            out.write("\tBoth reads in pair are mapped {0}\n".format(self.mm_count))
            out.write("\nStatistics of mapped pairs\n")
            out.write("\tSupplementary alignments {0}\n".format(self.supp_count))
            out.write("\tSecondary alignments {0}\n".format(self.sec_count))
            out.write("\tPrimary alignments {0}\n".format(self.prime_count))
            out.write("\nStatistics of primary alignments\n")
            out.write("\tWith mapq = 0 {0}\n".format(self.mapq0_count))
            out.write("\tWith mapq > 0 {0}\n".format(self.posq_count))
