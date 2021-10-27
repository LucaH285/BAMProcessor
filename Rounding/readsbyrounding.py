# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 16:33:56 2021

@author: Desktop
"""


class CountReadsByGrouping(Preprocessing):
    def ExportFxn2(self, Exprt, ExportLocation, ProcessedFrameLists, FileName):
        Export = Exprt.ExportFxn(ExportLocation, ProcessedFrameLists, FileName)
        return(Export)

    def VarLoads(self, Init):
        print("Initiated count by grouping...")
        Frames = BamFileLoads().VarLoads(Init)
        Genome = GenomeBinDivider().VarLoads(Init)
        ProcessedFrameList = []
        for DF in Frames:
            Sorted = DF.sort_values(by = ["Start", "End"])
            BinDict = {
                Bin: [0, []] for Bin in range(Init.GenomeStartRange, Init.GenomeEndRange) if Bin % Init.BinSize == 0
                }
            for Bin1, Bin2 in zip(Genome[:-1], Genome[1:]):
                for ReadStart, ReadEnd, Seq in zip(Sorted["Start"], Sorted["End"], Sorted["Sequence"]):
                    """
                    Condition 1, the read is fully within the Bin1, Bin2 iteration
                    Condition 2, the read starts in bin1, ends in bin2
                    """
                    if ((ReadStart >= Bin1) and (ReadStart < Bin2) and (ReadEnd > Bin1) and (ReadEnd <= Bin2)):
                        BinDict[Bin1][0] += 1
                        BinDict[Bin1][1].append(Seq)
                    elif ((ReadStart >= Bin1) and (ReadStart < Bin2) and (ReadEnd > Bin2)):
                        if (ReadEnd <= ((float(Init.Threshold) * Init.BinSize) + Bin2)):
                            BinDict[Bin1][0] += 1
                            BinDict[Bin1][1].append(Seq)
                        elif (ReadEnd >((float(Init.Threshold) * Init.BinSize) + Bin2)):
                            BinDict[Bin2][0] += 1
                            BinDict[Bin2][1].append(Seq)
            DataStructure = {
                "Genome Range Start Coordinate":[Bins for Bins in BinDict],
                "Reads per region":[BinDict[Reads][0] for Reads in BinDict],
                "Sequence":[BinDict[Seq][1] for Seq in BinDict]
                }
            ProcessedFrame = pd.DataFrame(DataStructure)
            ProcessedFrameList.append(ProcessedFrame)
        if Init.CreateAverageFrameFromFrames is True:
            AveragedFrame = [pd.concat([Frames for Frames in ProcessedFrameList]).
                             groupby(level=0).
                             mean(["Reads per region"])]
            if len(Init.ExportLocation) > 0:
                self.ExportFxn2(Exprt, Init.ExportLocation, AveragedFrame, FileName = "AveragedGroupedFrame")
            else:
                pass
            return(AveragedFrame)
        elif Init.CreateAverageFrameFromFrames is False:
            if len(Init.ExportLocation) > 0:
                self.ExportFxn2(Exprt, Init.ExportLocation, ProcessedFrameList, FileName = "GroupedFrame")
            else:
                pass
            return(ProcessedFrameList)