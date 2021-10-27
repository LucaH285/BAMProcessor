# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 16:33:50 2021

@author: Desktop
"""
class CountReadsDirectly(Preprocessing):
    def ExportFxn2(self, Exprt, ExportLocation, ProcessedFrameLists, FileName):
        Export = Exprt.ExportFxn(ExportLocation, ProcessedFrameLists, FileName)
        return(Export)

    def VarLoads(self, Init):
        Frames = BamFileLoads().VarLoads(Init)
        Genome = GenomeBinDivider().VarLoads(Init)
        ProcessedFrameList = []
        for DF in Frames:
            Sorted = DF.sort_values(by = ["Start", "End"])
            BinDict = {
                Bin: [0, []] for Bin in range(Init.GenomeStartRange, Init.GenomeEndRange) if Bin % Init.BinSize == 0
                }
            for Bin1, Bin2 in zip(Genome[:-1], Genome[1:]):
                Counter1 = 0
                Counter2 = 0
                for ReadStart, ReadEnd, Seq, ReadWidth in zip(Sorted["Start"], Sorted["End"],
                                                  Sorted["Sequence"], Sorted["Width"]):
                    if ((ReadStart >= Bin1) and (ReadStart < Bin2) and (ReadEnd > Bin1) and (ReadEnd <= Bin2)):
                        BinDict[Bin1][0] += 1
                        BinDict[Bin1][1].append(Seq)
                    elif ((ReadStart >= Bin1) and (ReadStart < Bin2) and (ReadEnd > Bin2)):
                        TempSeqList1 = []
                        TempSeqList2 = []
                        #Number of nucleotides in Bin1-Bin2 range
                        ReadInBlock = (Bin2 - ReadStart)
                        #Number of nucleotides outside the Bin1-Bin2 range
                        ReadOutBlock = (ReadEnd - Bin2)
                        #Fractional counts
                        Counter1 += ((ReadInBlock)/len(Seq))
                        Counter2 += ((ReadOutBlock)/len(Seq))
                        Condition = True
                        Count = 0
                        while(Condition):
                            if (Count < ReadInBlock):
                                TempSeqList1.append(Seq[Count])
                            elif ((Count >= ReadInBlock) and (Count <= len(Seq))):
                                try:
                                    TempSeqList2.append(Seq[Count])
                                except IndexError:
                                    pass
                            elif (Count > len(Seq)):
                                TempSeqList1 = "".join(TempSeqList1)
                                TempSeqList2 = "".join(TempSeqList2)
                                Condition = False
                            Count += 1
                        BinDict[Bin1][1].append(TempSeqList1)
                        BinDict[Bin2][1].append(TempSeqList2)
                BinDict[Bin1][0] += Counter1
                BinDict[Bin2][0] += Counter2
            DataStructure = {
                "Genome Range Start Coordinate":[Bins for Bins in BinDict],
                "Reads per region":[BinDict[Reads][0] for Reads in BinDict],
                "Sequence":[BinDict[Seq][1] for Seq in BinDict]
                }
            ProcessedFrame = pd.DataFrame(DataStructure)
            ProcessedFrameList.append(ProcessedFrame)
        if len(Init.ExportLocation) > 0:
            self.ExportFxn2(Exprt, Init.ExportLocation, ProcessedFrameList, FileName = "DirectCounts")
            return(ProcessedFrameList)
        else:
            return(ProcessedFrameList)