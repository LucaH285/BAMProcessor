# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 15:37:17 2021

@author: Desktop
"""

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from abc import ABC, abstractmethod



class AlignReads:
    def __init__(self, AlignmentStart, AlignmentEnd, ReferenceSequence, MainFrame):
            self.AlignStart = AlignmentStart
            self.AlignEnd = AlignmentEnd
            self.RefSeq = ReferenceSequence
            self.InputFile = MainFrame
         
class Alignments(ABC):
    @abstractmethod
    def ReadAlign(self, Align):
        pass
            
class ReadAlignment(Alignments):
    def ReadAlign(self, Align):
        InputFile = Align.InputFile
        Start = Align.AlignStart
        End = Align.AlignEnd
        RefSeq = Align.RefSeq
        for Frames, rng in zip(InputFile, range(len(InputFile))):
            ReadSeq = [Frames["Sequence"][Ind]
                       for Ind in Frames.index.values
                       if Frames["Start"][Ind] >= Start
                       and Frames["End"][Ind] <= End]
            AlignmentsList = []
            #Align reads to each other
            if len(RefSeq) == 0:
                print("No refseq")
                for Seqs1, Seqs2 in zip(ReadSeq[:-1], ReadSeq[1:]):
                    alignments = pairwise2.align.globalms(Seqs1, Seqs2, 2, -2, -1, -.1)
                    AlignmentsList.append(format_alignment(*alignments[0]))
                    print(format_alignment(*alignments[0]))
            #Align reads to a reference sequence
            elif len(RefSeq) != 0:
                for Seqs in ReadSeq:
                    alignments = pairwise2.align.globalms(RefSeq, Seqs, 2, -2, -1, -.1)
                    AlignmentsList.append(format_alignment(*alignments[0]))
                    print(format_alignment(*alignments[0]))
            print("For frame {0}, alignments, {1}".format(rng, AlignmentsList[0:]))

