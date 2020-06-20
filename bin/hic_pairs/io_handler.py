import logging
from typing import List, Tuple

import pysam

from gfa_graph.asm_graph import AsmGraph

logger = logging.getLogger()

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

        # total number of reads
        self.read_counter = 0
        # total number of alignment records
        self.alignment_counter = 0
        # total number of reads that does not have primary alignments == should be zero
        self.wo_primary_alignments = 0
        # n = non-mapped, m = mapped, each letter is read from pair
        self.nn_counter, self.nm_counter, self.mm_counter = 0, 0, 0
        # supp = supplementary, sec = secondary, prime = prim
        self.supp_counter, self.sec_counter, self.prime_counter = 0, 0, 0
        # mapq = 0 is interesting since it may have perfect alignments but with XA tag
        self.mq00_counter, self.mq01_counter, self.mq11_counter = 0, 0, 0
        # posq can be further filtered
        self.posq_counter = 0
        # type of hic pairs
        self.uu_pair_counter, self.um_pair_counter, self.mm_pair_counter = 0, 0, 0

        self.very_bad = 0
        self.uniquely_mapped_reads = dict()

    # def parse_primary_mapped_reads(self, info):
    #     if 0 != info.query_alignment_start:
    #         if info.flag & 0x0010 == 0:
    #             strand = '+'
    #             print(info.to_string())
    #             print("Positive {0} {1} {2}".format(info.reference_start, info.reference_end, info.cigarstring))
    #         else:
    #             strand = '-'
    #             print(info.to_string())
    #             print("Negative {0} {1} {2}".format(info.reference_start, info.reference_end, info.cigarstring))

    def _partition_read_alignments(self, records: List) -> Tuple[List, List, List]:
        primary, secondary, supplementary = [], [], []

        for record in records:
            if record.flag & 0x0100:  # info.is_secondary
                secondary.append(record)
            elif record.flag & 0x0800:  # info.is_supplementary
                supplementary.append(record)
            else:
                primary.append(record)

        if len(primary) != 1:
            return [], [], []

        return primary, secondary, supplementary

    def _parse_paired_read_alignments(self, recs1: List, recs2: List) -> None:
        if len(recs1) == 0 or len(recs2) == 0:
            return

        prim_rs1, sec_rs1, supp_rs1 = self._partition_read_alignments(recs1)
        prim_rs2, sec_rs2, supp_rs2 = self._partition_read_alignments(recs2)

        if len(prim_rs1) != 1 or len(prim_rs2) != 1:
            self.wo_primary_alignments += 1
            logger.error("It must be one and only one primary alignment for both reads")
            return

        prim_r1, prim_r2 = prim_rs1[0], prim_rs2[0]

        if prim_r1.flag & 0x0004 and prim_r2.flag & 0x0004:
            # both reads are unmapped - we are not interested
            self.nn_counter += 1
            return
        elif bool(prim_r1.flag & 0x0004) != bool(prim_r2.flag & 0x0004):
            # one of two reads is unmapped - we are not interested
            self.nm_counter += 1
            return

        self.mm_counter += 1

        if not (len(sec_rs1) or len(sec_rs2) or len(supp_rs1) or len(supp_rs2)):
            self.prime_counter += 1

            if prim_r1.has_tag('SA') or prim_r2.has_tag('SA'):
                logger.error("There is a tag SA, but there is no records for that.")
                return

            if prim_r1.mapping_quality == 0 and prim_r2.mapping_quality == 0:
                # MM type
                self.mq00_counter += 1

                if not prim_r1.has_tag('XA') and not prim_r2.has_tag('XA'):
                    self.very_bad += 1
                    # logger.debug("There is no XA tags and mapq = 0, there is nothing we can do it here.")
                    return

            elif (prim_r1.mapping_quality == 0) != (prim_r2.mapping_quality == 0):
                # One of alignments unique, another one is not
                # assert prim_r1.has_tag('XA') or prim_r2.has_tag('XA')
                self.mq01_counter += 1
            else:
                # we have unique mapping - hooray!
                self.mq11_counter += 1

        if len(sec_rs1) or len(sec_rs2):
            self.sec_counter += 1

        if len(supp_rs1) or len(supp_rs2):
            self.supp_counter += 1

    def parse(self) -> None:
        """
        The main assumption is the bam file sorted by read name
        """
        recs1, recs2 = [], []

        def push_alignments(rec) -> None:
            if rec.flag & 0x0040:  # rec.is_read1
                recs1.append(rec)
            elif rec.flag & 0x0080:  # rec.is_read2
                recs2.append(rec)
            else:
                assert False

        prev_read_id = None
        bam_stream = pysam.AlignmentFile(self.path_to_bam, 'rb')
        for record in bam_stream.fetch(until_eof=True):
            self.alignment_counter += 1

            read_id = record.query_name

            if not prev_read_id:
                prev_read_id = read_id
                push_alignments(record)
                continue

            if read_id != prev_read_id:
                # logger.debug('New group of alignments for {0} is detected.'.format(read_id))
                self.read_counter += 1

                self._parse_paired_read_alignments(recs1, recs2)

                recs1.clear()
                recs2.clear()
                prev_read_id = read_id

            push_alignments(record)

            if self.read_counter % 100000 == 0:
                pass
                # logger.info("{0} reads were processed.".format(self.read_counter))

        if len(recs1) != 0 or len(recs2) != 0:
            # logger.info("Dump the latest group alignment")
            self.read_counter += 1
            self._parse_paired_read_alignments(recs1, recs2)

    def write_statistics_info_into_file(self, out_file: str) -> None:
        with open(out_file, 'w') as out:
            out.write("Total number of reads = {0}\n".format(self.read_counter))
            out.write("Total number of alignment records = {0}\n".format(self.alignment_counter))
            out.write("\tBoth reads in pair are unmapped {0}\n".format(self.nn_counter))
            out.write("\tOne of the reads in pair is unmapped {0}\n".format(self.nm_counter))
            out.write("\tBoth reads in pair are mapped {0}\n".format(self.mm_counter))
            out.write("\tReads without any or more than one "
                      "primary alignments {0}\n".format(self.wo_primary_alignments))
            out.write("\nStatistics of alignment types\n")
            out.write("\tSupplementary alignments {0}\n".format(self.supp_counter))
            out.write("\tSecondary alignments {0}\n".format(self.sec_counter))
            out.write("\tPrimary alignments {0}\n".format(self.prime_counter))
            out.write("\nStatistics of primary alignments\n")
            out.write("\tBoth reads have mapq = 0 {0}\n".format(self.mq00_counter))
            out.write("\tBoth reads have mapq = 0 and no XA tags {0}\n".format(self.very_bad))
            out.write("\tOne of the reads have mapq = 0 {0}\n".format(self.mq01_counter))
            out.write("\tBoth reads have mapq > 0 {0}\n".format(self.mq11_counter))
            # out.write("\tWith mapq > 0 {0}\n".format(self.posq_counter))
            out.write("\n\nStatistics of Hi-C pairs\n")
            out.write("\tUU type {0}\n".format(self.uu_pair_counter))
            out.write("\tUM type {0}\n".format(self.um_pair_counter))
            out.write("\tMM type {0}\n".format(self.mm_pair_counter))


# print(info)
# qname = info.query_name  # info.
# rname = info.reference_name
# pos = info.reference_start
# cigar = info.cigarstring
# mrnm = info.next_reference_name
# mpos = info.next_reference_start
# tags = 0
# if 0 != info.query_alignment_start:
#     print(mrnm + " " + cigar)
#     print(str(pos) + " " + str(info.query_alignment_start))
#     # self.uniquely_mapped_reads

# # self.mm_counter += 1
# if info.flag & 0x0800:  # info.is_supplementary:
#     # TODO small overlaps because of repeats or chimeric stuff
#     pass
# elif info.flag & 0x0100:  # info.is_secondary:
#     # TODO may cause by overlaps in unitigs or complete overlaps
#     pass
# else:
#     self.prime_counter += 1
#     if info.mapping_quality == 0:
#         self.mapq0_counter += 1
#     else:
#         self.posq_counter += 1
#         # self.parse_primary_mapped_reads(info)

# if len(recs1) == 1 and len(recs2) == 1:
#     if recs1[0].flag & 0x0004 and recs2[0].flag & 0x0004:
#         # both reads are unmapped - we are not interested
#         self.nn_counter += 1
#     elif bool(recs1[0].flag & 0x0004) != bool(recs2[0].flag & 0x0004):
#         # one of two reads is unmapped - we are not interested
#         self.nm_counter += 1
#     elif recs1[0].mapping_quality > 0 and recs2[0].mapping_quality > 0:
#         # assert recs1[0].flag & 0x0800 == 0 and recs2[0].flag & 0x0100 == 0
#         # logger.debug("We have uniquely mapped pair")
#         self.uu_pair_counter += 1
#     elif bool(recs1[0].mapping_quality) != bool(recs2[0].mapping_quality):
#         # assert recs1[0].flag & 0x0800 == 0 and recs2[0].flag & 0x0100 == 0
#         self.um_pair_counter += 1
#     else:
#         assert recs1[0].mapping_quality == 0 and recs2[0].mapping_quality == 0
#         if recs2[0].flag & 0x0100 != 0:
#             print(recs2[0].tostring())
#         # assert recs1[0].flag & 0x0800 == 0 and recs2[0].flag & 0x0100 == 0
#         self.mm_pair_counter += 1
# else:
#     pass
# def _parse_read_alignments(self, records: List):
#     for info in records:
#         self.alignment_counter += 1
#         if info.flag & 0x0004 and info.flag & 0x0008:
#             # both reads are unmapped - we are not interested
#             self.nn_counter += 1
#             # pass
#         elif bool(info.flag & 0x0004) != bool(info.flag & 0x0008):
#             # one of two reads is unmapped - we are not interested
#             self.nm_counter += 1
#             # pass
#         else:
#             self.mm_counter += 1
#             if info.flag & 0x0800:  # info.is_supplementary:
#                 # TODO small overlaps because of repeats or chimeric stuff
#                 self.supp_counter += 1
#                 # pass
#             elif info.flag & 0x0100:  # info.is_secondary:
#                 # TODO may cause by overlaps in unitigs or complete overlaps
#                 self.sec_counter += 1
#                 # pass
#             else:
#                 self.prime_counter += 1
#                 if info.mapping_quality == 0:
#                     self.mapq0_counter += 1
#                 else:
#                     self.posq_counter += 1
#                     # self.parse_primary_mapped_reads(info)
