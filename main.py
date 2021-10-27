# -*- coding: utf-8 -*-
"""
Created on Sun Jul 25 16:24:29 2021

@author: Luca Hategan
@Description: ATAC-seq signal filtering, and PcCRE detection
@Version: 2.0
"""
#######################
#Library imports
#######################
import pandas as pd
import numpy as np

import bamnostic as bs
import os
from abc import ABC, abstractmethod
import itertools
import argparse
import warnings
import distutils
from scipy.stats import nbinom
#Only for debugging:
import time

###################
#A method constant variable
#to be used only when wanting
#to test/run files through the 
#IDE instead of the CLI.
#
#Can be useful if you want to 
#test your files. 
#Default is to use the CLI
###################
VARLOAD_METHOD = "useIDE"

CWD = os.getcwd()

class InitializeVars:
    def __init__(self, List, Start, End,
                 Chromosome, BinSize, Thresh,
                 Strnd, Export, CCREs, AccountForOpenSpread,
                 Bool, AverageFrames):
        self.Bam_file_list = List
        self.GenomeStartRange = Start
        self.GenomeEndRange = End
        self.Chr = Chromosome
        self.BinSize = BinSize
        self.Threshold = Thresh
        self.Strand = Strnd
        self.ExportLocation = Export
        self.KnownCCRE = CCREs
        self.RangeAccount = AccountForOpenSpread
        self.CTCFSiteAsFirstCCRE = Bool
        self.CreateAverageFrameFromFrames = AverageFrames
        self.InputFrames = []
        self.BinnedGenome = []
            
class Sorters:
    def GeneralSorting(self, List):
        Condition = True
        while(Condition):
            Condition = False
            for Ind in range(len(List) - 1):
                if List[Ind] > List[Ind + 1]:
                    List[Ind], List[Ind + 1] = List[Ind + 1], List[Ind]
                    Condition = True
        return(List)

class ExportClass:
    def ExportFxn(self, ExportLocation, ProcessedFrameList, FileName):
        if len(ExportLocation) < len(ProcessedFrameList):
            SizeList = [len(ProcessedFrameList), len(ExportLocation)]
            SizeDifference = max(SizeList) - min(SizeList)
            #Export Location + the final element in the Export Location List repeated for
            #n times, where n is the difference between the ProcessedFrame inputs and the
            #Export location.
            NewExportList = ExportLocation + [ExportLocation[len(ExportLocation) - 1] for _ in range(SizeDifference)]
            for Exports, Frames, rng in zip(NewExportList, ProcessedFrameList, range(len(NewExportList))):
                ExportName = "{0}_{1}.csv".format(FileName, rng)
                ExportStmnt = Exports + ExportName
                Frames.to_csv(ExportStmnt)
        elif len(ExportLocation) == len(ProcessedFrameList):
            for Exports, Frames, rng in zip(ExportLocation, ProcessedFrameList, range(len(ExportLocation))):
                ExportName = "{0}_{1}.csv".format(FileName, rng)
                ExportStmnt = Exports + ExportName
                Frames.to_csv(ExportStmnt)

class Statistics:
    def TwoSigmaLOD(self, InputFrame, KnownCCRE):
        ReadsPerRegionStdDev = np.std(InputFrame["Reads per region"])
        IndexVals = []
        StdDev_Structure = {
            "Std Dev":[ReadsPerRegionStdDev],
            "2-Sigma":[2*ReadsPerRegionStdDev]
            }
        Frame = pd.DataFrame(StdDev_Structure)
        for Counts, CountRange in zip(InputFrame["Reads per region"], InputFrame.index.values):
            if Counts > (2 * ReadsPerRegionStdDev):
                pass
            else:
                #Here, cannot use this way of retrieving the index values, since not every
                #Reads per region entry is unique. Ie: there are many index values that have
                #an associated Reads per region = 2
                #print(InputFrame.index[InputFrame["Reads per region"] == Counts], Counts)
                IndexVals.append(CountRange)
        PcCREs = InputFrame.drop([Ind for Ind in IndexVals], axis = 0).reset_index(inplace = False, drop = True)
        ClosestCCRE = []
        Distance = []
        #The issue here is that both PCCRE frames from two seperate bam files are being processed at the same time.
        for PCCRE in PcCREs["Genome Range Start Coordinate"]:
            """
            Notice that while in the loop, KCCRE stays constant untill the next loop point,
            really only need to discriminate between the X and Y which represent start and end
            coordinates of the known CCRE.

            Admitidely this is a very inefficient way to do this (probably)

            Initiates KCCRE 1...
            """
            KCCRE_X = []
            KCCRE_Y = []
            KCCRE_Start = []
            KCCRE_End = []
            for KCCRE in KnownCCRE:
                """
                loops through all PCCREs for each KCCRE
                i.e.: KCCRE1:PCCRE1, KCCRE1:PCCRE2, KCCRE1:PCCRE3...
                """
                KCCRE_X.append(abs(PCCRE - KCCRE[0]))
                KCCRE_Y.append(abs(PCCRE - KCCRE[1]))
                KCCRE_Start.append(KCCRE[0])
                KCCRE_End.append(KCCRE[1])
            #These minimums match the index value of the KCCRE coordinates
            #So I retrieve the correct KCCRE by using the index of the minimum value
            #in the KCCRE_X list, which houses the difference between PCCRE and KCCRE
            MinimumsList = [min(KCCRE_X), min(KCCRE_Y)]
            GetIndex_X = KCCRE_X.index(MinimumsList[0])
            GetIndex_Y = KCCRE_Y.index(MinimumsList[1])
            Distance.append(min(MinimumsList))
            if min(MinimumsList) == MinimumsList[0]:
                #The index is then used to retrieve the value of the KCCRE coordinate from the KCCRE_Start list
                ClosestCCRE.append((KCCRE_Start[GetIndex_X], "kccre_start_coordinate"))
            elif min(MinimumsList) == MinimumsList[1]:
                ClosestCCRE.append((KCCRE_Start[GetIndex_Y], "kccre_end_coordinate"))
        ClosestCCREDataStructure = {
            "closestCCRE":ClosestCCRE,
            "distance":Distance
            }
        ClosestCCREFrame = pd.DataFrame(ClosestCCREDataStructure)
        Concat = pd.concat([PcCREs, ClosestCCREFrame, Frame], axis = 1)
        Concat = Concat.fillna("")
        return(Concat)

    def negativeBinomial(self, InputFrame):
        Mu = np.average(InputFrame["Reads per region"])
        Var = np.var(InputFrame["Reads per region"])
        Prob = Mu/Var
        N = ((Mu ** 2)/(Var - Mu))
        MinK = min(InputFrame["Reads per region"])
        MaxK = max(InputFrame["Reads per region"])
        ProbabilitVector = [nbinom.pmf(range(MinK, MaxK), n = N, p = Prob)]
        return("Not yet finished")
    
    def normalityModel(self, InputFrame):
        """
        Call with direct 
        """
        print(InputFrame)
        breakpoint()

class Preprocessing(ABC):
    @abstractmethod
    def VarLoads(self, Init):
        pass

class InheritableSorting(ABC):
    @abstractmethod
    def Sorting(self, Sort):
        pass



class Exports(ABC):
    @abstractmethod
    def ExportFxn2(self, Exprt):
        pass

class Stats(ABC):
    @abstractmethod
    def StatisticalFunctions(self, Stats):
        pass

class BamFileLoads(Preprocessing):
    def VarLoads(self, Init):
        print("Processing...")
        global FrameList
        FrameList = []
        DictionaryList = []
        #First create the index files, if needed
        for Files in range(len(Init.Bam_file_list)):
            # global CWD
            # CWD = os.getcwd()
            Split = Init.Bam_file_list[Files].split('\\')
            os.chdir(Split[0])
            Bam = bs.AlignmentFile(Init.Bam_file_list[Files], "rb")
            try:
                BamDict = {"read_lis": [], "read_info": [], "mapq_lis": [], "seq_lis": [], "width_lis": []}
                for Reads in Bam.fetch(Init.Chr, Init.GenomeStartRange, Init.GenomeEndRange):
                    BamDict["read_lis"].append((Reads.reference_start, Reads.reference_end))
                    BamDict["read_info"].append(Reads.is_reverse)
                    BamDict["mapq_lis"].append(Reads.mapq)
                    BamDict["seq_lis"].append(Reads.seq)
                    BamDict["width_lis"].append(Reads.reference_length)
                DictionaryList.append(BamDict)
            except ValueError:
                print("No bai file found for {}, creating bai index file...".format(Init.Bam_file_list[Files]))
                try:
                    from rpy2 import robjects
                    from rpy2.robjects.packages import importr
                    import rpy2.rinterface
                    base = importr('base')
                    utils = importr('utils')
                    #Adjust Path to be compatible with R
                    Path = Init.Bam_file_list[Files].split(sep="\\")
                    PathJoined = "/".join(Path)
                    Glb = robjects.globalenv['Object'] = PathJoined
                    Lib = robjects.r('library(Rsamtools)')
                    Fxn = robjects.r('indexBam(Object)')
                except:
                    raise(TypeError("""
                        Please either install rpy2 from anaconda, and if you've done so and this error appears, or JIT failure error, 
                        please add rpy2 to your %PATH%, see: https://jianghaochu.github.io/how-to-install-rpy2-in-windows-10.html"
                        """))
        for Dicts in DictionaryList:
            Structure = {
                "Mapq":Dicts["mapq_lis"],
                "Sequence":Dicts["seq_lis"],
                "Start":[i[0] for i in Dicts["read_lis"]],
                "End":[j[1] for j in Dicts["read_lis"]],
                "Width":Dicts["width_lis"],
                "Strand Order":Dicts["read_info"]
                }
            BamAsFrame = pd.DataFrame(Structure)
            for rows, strand in zip(BamAsFrame.index, BamAsFrame["Strand Order"]):
                if Init.Strand == "+":
                    if str(strand) != "False":
                        BamAsFrame.drop(rows, axis = 0, inplace = True)
                elif Init.Strand == "-":
                    if str(strand) != "True":
                        BamAsFrame.drop(rows, axis = 0, inplace = True)
                elif Init.Strand == "*":
                    pass
            FrameList.append(BamAsFrame)
        return(FrameList)

class GenomeBinDivider(Preprocessing):
    def VarLoads(self, Init):
        Counter = Init.GenomeStartRange
        EndRange = Init.GenomeEndRange
        BinnedGenome = []
        while Counter < EndRange:
            if Counter == Init.GenomeStartRange:
                BinnedGenome.append(Counter)
                Counter += Init.BinSize
            else:
                BinnedGenome.append(Counter)
                Counter += Init.BinSize
        return(BinnedGenome)




        
class GetReadCountsPerCCRE(Preprocessing, InheritableSorting, Exports):
    def Sorting(self, Sort, List1):
        Sorting = Sort.GeneralSorting(List1)
        return(Sorting)

    def ExportFxn2(self, Exprt, ExportLocation, ProcessedFrameLists, FileName):
        Export = Exprt.ExportFxn(ExportLocation, ProcessedFrameLists, FileName)
        return(Export)

    def VarLoads(self, Init, InputFiles):
        Frames = FrameList
        KnownCCREs1 = Init.KnownCCRE
        SortedX = self.Sorting(Sort, [xCoor[0] for xCoor in KnownCCREs1])
        SortedY = self.Sorting(Sort, [yCoor[1] for yCoor in KnownCCREs1])
        KnownCCREs = [(X, Y) for X, Y in zip(SortedX, SortedY)]
        ProcessedFrameList = []
        for frames in Frames:
            CCREDict = {
                CCRECoord: 0 for CCRECoord in KnownCCREs
                }
            for CCREs, rng in zip(KnownCCREs, range(len(KnownCCREs) - 1)):
                CCREStart = CCREs[0]
                CCREEnd = CCREs[1]
                for Start, End in zip(frames["Start"], frames["End"]):
                    if ((Start >= CCREStart) and (End <= CCREEnd)):
                        CCREDict[(CCREStart, CCREEnd)] += 1
                    elif ((Start >= CCREStart) and (End > CCREEnd)):
                        Width = (End - Start)
                        if (End <= (((Width) * 0.5) + CCREEnd)):
                            CCREDict[(CCREStart, CCREEnd)] += 1
            DataStructure = {
                "CCRE Coordinates":[Coors for Coors in CCREDict],
                "Reads per region":[CCREDict[Reads] for Reads in CCREDict]
                }
            AverageStdDevStructure = {
                "Average":[np.average([CCREDict[Reads] for Reads in CCREDict])],
                "StdDev":[np.std([CCREDict[Reads] for Reads in CCREDict])]
                }
            CCRECountFrame = pd.DataFrame(DataStructure)
            CCREAvgStdDevTempFrame = pd.DataFrame(AverageStdDevStructure)
            ExportFrame = pd.concat([CCRECountFrame, CCREAvgStdDevTempFrame], axis = 1)
            ExportFrame = ExportFrame.fillna(value = "")
            ProcessedFrameList.append(ExportFrame)
        if len(Init.ExportLocation) > 0:
            self.ExportFxn2(Exprt, ExportLocation=Init.ExportLocation, ProcessedFrameLists=ProcessedFrameList, FileName = "CountsPerCCRE")
        else:
            pass
        return(ProcessedFrameList)

class RemoveKnownCCRE(Preprocessing, Exports, InheritableSorting):
    def ExportFxn2(self, Exprt, ExportLocation, ProcessedFrameLists, FileName):
        Export = Exprt.ExportFxn(ExportLocation, ProcessedFrameLists, FileName)
        return(Export)

    def Sorting(self, Sort, List1):
        Sorting = Sort.GeneralSorting(List1)
        return(Sorting)

    def VarLoads(self, Init, InputFiles):
        print("Active")
        KnownCCREs1 = Init.KnownCCRE
        SortedX = self.Sorting(Sort, [xCoor[0] for xCoor in KnownCCREs1])
        SortedY = self.Sorting(Sort, [yCoor[1] for yCoor in KnownCCREs1])
        KnownCCREs = [(X, Y) for X, Y in zip(SortedX, SortedY)]
        IndexList = []
        FrameList = []
        if InputFiles is not None:
            for InFrames in InputFiles:
                for cCREs, rng in itertools.product(range(len(KnownCCREs)), InFrames.index.values):
                    StartCoor = InFrames["Genome Range Start Coordinate"][rng]
                    #This difference function calculates the distance
                    #Between the known cCRE and the Bin coordinates
                    DiffX = abs(KnownCCREs[cCREs][0] - StartCoor)
                    DiffY = abs(KnownCCREs[cCREs][1] - StartCoor)
                    # print(KnownCCREs[cCREs][0], KnownCCREs[cCREs][1], StartCoor)
                    # time.sleep(1)
                    if Init.CTCFSiteAsFirstCCRE is False:
                        if Init.RangeAccount > 0:
                            if ((DiffX <= Init.RangeAccount) or (DiffY <= Init.RangeAccount)):
                                IndexList.append(rng)
                        elif Init.RangeAccount == 0:
                            if ((StartCoor >= KnownCCREs[cCREs][0]) and (StartCoor <= KnownCCREs[cCREs][1])):
                                IndexList.append(rng)
                    elif Init.CTCFSiteAsFirstCCRE is True:
                        if cCREs == 0:
                            if ((StartCoor >= KnownCCREs[cCREs][0]) and (StartCoor <= KnownCCREs[cCREs][1])):
                                IndexList.append(rng)
                        elif cCREs > 0:
                            if Init.RangeAccount > 0:
                                if ((DiffX <= Init.RangeAccount) or (DiffY <= Init.RangeAccount)):
                                    IndexList.append(rng)
                            elif Init.RangeAccount == 0:
                                if ((StartCoor >= KnownCCREs[cCREs][0]) and (StartCoor <= KnownCCREs[cCREs][1])):
                                    IndexList.append(rng)
            for Frames in InputFiles:
                if "Sequence" in Frames.columns:
                    Frames = Frames.drop(columns = ["Sequence"])
                    Frames = Frames.drop([i for i in IndexList], axis = 0)
                else:
                    Frames = Frames.drop([i for i in IndexList], axis = 0)
                FrameList.append(Frames)
            if len(Init.ExportLocation) > 0:
                self.ExportFxn2(Exprt, Init.ExportLocation, ProcessedFrameLists=FrameList, FileName="RemovedKnownCCREs")
            else:
                pass
        else:
            pass
        return(FrameList)

class StatisticalProcessing(Stats, Exports):
    def ExportFxn2(self, Exprt, ExportLocation, ProcessedFrameLists, FileName):
        Export = Exprt.ExportFxn(ExportLocation, ProcessedFrameLists, FileName)
        return(Export)

    def StatisticalFunctions(self, Stats, Inputs):
        PcCREFrames = []
        for Frames in Inputs:
            Processing = Stats.TwoSigmaLOD(Frames, KnownCCRE=Init.KnownCCRE)
            PcCREFrames.append(Processing)
        if len(Init.ExportLocation) > 0:
            self.ExportFxn2(Exprt, Init.ExportLocation, ProcessedFrameLists=PcCREFrames, FileName="PotentialCCREs")
        else:
            pass
        return(PcCREFrames)



class ResetCWDToOld:
    """
    This class, when called, resets the working directory to the location directory.
    I cannot think of any way to circumvent this since an error is raised if the
    imports are not in the same directory/drive as the inital CWD.

    So the CWD is changed in the bam loads class and must be set back to its original
    value when the programs ends.
    """
    def ResetFxn(self, InputCWD):
        return(os.chdir(InputCWD))

class argumentParser:
    def useParser(self):
        def CSVConverter(Input, Seperation=",", convType = str):
            #Should raise an error here
            Converted = []
            for Vals in Input.split(Seperation):
                Converted.append(convType(Vals))
            return(Converted)
            
        def CSVConvertIntegers(Input, Seperation=",", convType = int):
            Converted = []
            if (Input.endswith(".csv") or Input.endswith(".txt")):
                Fxn = [str(Vals) for Vals in Input.split(",")]
                return(Fxn)
            else:
                for Vals in Input.split(Seperation):
                    Converted.append(convType(Vals))
            return(Converted)
        
        parser = argparse.ArgumentParser(description="ATAC-Seq processor, Version 2.0")
        parser.add_argument("-l", "--list", help="add Bam file paths",
                            required=True, type=CSVConverter, dest="files")
        parser.add_argument("-S", "--StartCoor", type=int, help="Start coordinate of the desired query region",
                            required=True, dest="StartCoordinate")
        parser.add_argument("-E", "--EndCoor", type=int, help="End coordinate of the desired query region", required=True, dest="EndCoordinate")
        parser.add_argument("-Chr", "--Chromosome", type=str, help="Chromosome of the desired query region, input as: 'chr2'",
                            required=True, dest="Chromosome")
        parser.add_argument("-Strnd", "--Strand", type=str, help="The strand of interest, +, - or *", required = True, dest = "Strand")
        parser.add_argument("-BS", "--BinSize", metavar="N", type=int, help="Size of the bins in nucleotides. This will be how the query regions is split",
                            required=False, dest="BinSize")
        parser.add_argument("-CCRES", "--CCREStartCoord", type=CSVConvertIntegers, help="Start Coordinates of the known CCRE's from ENCODE", dest="CCREStart")
        parser.add_argument("-CCREE", "--CCREEndCoord", type=CSVConvertIntegers, help="End Coordinates of the known CCRE's from ENCODE", dest="CCREEnd")
        parser.add_argument("-Exprt", "--ExportLocation", action="append", help="add export location paths + tissue and PND identifier i.e.: C:\Test\MidbrainPND0",
                            required=False, metavar="FILE", dest="ExportFile")
        parser.add_argument("-Sprd", "--Spread", metavar="N", type=int, help="Value accounting for the possible spread of open chromatin surrounding a CCRE, in nt",
                            required=False, dest="SpreadAccount", default=0)
        parser.add_argument("-Thresh", "--Threshold", metavar="N", type=float,
                            help="Value for the the CounReadsByGrouping function, assigns how the fraction of the read needed to overlap with an adjacent bin in order to be a part of it",
                            required=True, dest="Threshold", default=0.33)
        parser.add_argument("-CTCF", "--FirstCCREasCTCFSite",
                            help="Processes the first CCRE as a CTCF binding site. May be important for genes with small 5' intergenic regions. Ignores the --Spread value for the first CCRE",
                            required=False, dest="BoolValCTCF", type=lambda x:bool(distutils.util.strtobool(x)), default=False)
        parser.add_argument("-Avg", "--AverageFrames",
                            help="Average the input files", type=lambda x:bool(distutils.util.strtobool(x)),
                            dest="AvgFrames", required=False, default=False)
        parser.add_argument("-Counts", "--GroupingFunction", type=str, help="Select the grouping function to use, either 'direct' or 'rounding'", required=True,
                            dest = "CountSelection")
        parser.add_argument("-A", "--Align", required=True,
                            help="string to perform alignment on reads: 'y' or 'n'", dest = "PerformAlign", type=str)
        parser.add_argument("-SA", "--StartAlign", metavar="N", type=int, help="Start coordinate of the alignment region",
                            required=False, dest="StartAlign")
        parser.add_argument("-EA", "--EndAlign", metavar="N", type=int, help="End coordinate of the alignment region",
                            required=False, dest="EndAlign")
        parser.add_argument("-Ref", "--ReferenceSequence", type=str, help="A string of the reference sequence to align reads to, i.e.: 'ATTCGACGGGTA'",
                            required=False, dest="RefSeq", default = [])
        args = parser.parse_args()
        return(args)
    
    def useIDE(self):
        print("Run Function Initiated")
        pass
        
if __name__ == '__main__':
    arg = argumentParser()
    if VARLOAD_METHOD == "useCLI":
        args = arg.useParser()
        for files in args.files:
            ext = os.path.splitext(files)[-1].lower()
            if str(ext) != ".bam":
                raise(ValueError("Bam file not entered, please check your inputs. The only NGS data structures currently supported are files ending with .bam"))
            if not os.path.exists(files):
                raise(FileNotFoundError("{0} Does not exist, please check your input".format(files)))
        
        BamFiles = args.files
        StartCoordinate = args.StartCoordinate
        EndCoordinate = args.EndCoordinate
        Chromosome = args.Chromosome
        if args.BinSize is None:
            BinSize = 100
            warnings.warn("No BinSize value entered, set to default: 100nt")
        else:
            BinSize = args.BinSize
    
        if (args.CCREStart is None) or (args.CCREEnd is None) :
            warnings.warn("Either No CCRE's entered, or missing values for CCRE Start/End Coordinates")
            KnownCCRE = [(0, 0)]
        elif ((len(args.CCREStart) != len(args.CCREEnd))):
            raise(ValueError("Coordinate mismatch in the inputted CCREs, check inputs. Make sure that the start and end coordiantes are entered correctly"))
        else:
            KnownCCRE = [(x, y) for x, y in zip(args.CCREStart, args.CCREEnd)]    
            
        if args.ExportFile is not None:
            ExportFile = args.ExportFile
        else:
            warnings.warn("No Export file(s) inputted")
            ExportFile = []
        Spread = args.SpreadAccount
        Threshold = args.Threshold
        if args.BoolValCTCF is None:
            warnings.warn("No Bool val for CTCF entered, default is False")
            Boolean = False
        else:
            Boolean = bool(args.BoolValCTCF)
        if args.AvgFrames is None:
            warnings.warn("No Bool val for Average Frames, default is False")
            AverageFrames = False
        else:
            AverageFrames = bool(args.AvgFrames)
        Strand = args.Strand
        Init = InitializeVars(
            List = BamFiles, Start = StartCoordinate, Strnd=Strand, Bool = Boolean,
            End = EndCoordinate, Chromosome = Chromosome, BinSize = BinSize, AverageFrames = AverageFrames,
            Thresh = Threshold, Export = ExportFile, CCREs = KnownCCRE, AccountForOpenSpread = Spread
            )
        
        Sort = Sorters()
        Exprt = ExportClass()
        Stats = Statistics()
        
        if str(args.CountSelection) != "direct" and str(args.CountSelection) != "rounding":
            raise(NameError("Please chose either 'direct' or 'rounding'. Ensure all characters are lower case"))
        else:
            global Process
            if str(args.CountSelection) == "direct":
                Process = CountReadsDirectly().VarLoads(Init)
            elif str(args.CountSelection) == "rounding":
                # global Process
                Process = CountReadsByGrouping().VarLoads(Init)
        ReadCountsCCRE = GetReadCountsPerCCRE()
        ReadCountsCCRE.VarLoads(Init, InputFiles = Process)
        
        Sts = StatisticalProcessing()
        Sts.StatisticalFunctions(Stats, Inputs = RemoveKnownCCRE().VarLoads(Init, InputFiles = Process))
        
        if args.PerformAlign == "y":
            from Test.Align import aligner
            Start = args.StartAlign
            End = args.EndAlign
            RefSeq = args.RefSeq
            Align = aligner.AlignReads(
                AlignmentStart=Start, 
                AlignmentEnd=End, 
                ReferenceSequence=RefSeq,
                MainFrame = FrameList
                )
            AlignFxn = aligner.ReadAlignment()
            AlignFxn.ReadAlign(Align)
            
        elif args.PerformAlign == "n":
            pass
        else:
            raise(TypeError("Please enter either 'y' for yes or 'n' in the -A option for alignments"))
    
        if "CWD" in globals():
            ResetCWDToOld().ResetFxn(InputCWD=CWD)
            
    elif VARLOAD_METHOD == "useIDE":
        args = arg.useIDE()
        ###################
        #Add in your inputs here
        #Default values are already set, empty strings or 0 values are vars to be populated
        ###################
        #Set CCREs here, in any order
        CCREList = [22618059,22620124,22620525,22620925,22621289,22621644,22621971,22622347,22623005,22623238,22623473,
                             22623833,22624142,22624636,22625172,22625544,22627266,22627922,22634813,22646161,22671008,22671524,
                             22618244,22620282,22620867,22621232,22621628,22621942,22622228,22622624,22623225,22623466,22623688,22624111,
                             22624360,22624897,22625505,22625892,22627514,22628230,22635078,22646505,22671251,22671685]
        
        CCREList.sort()
        CCRE_StartCoords = [CCREList[x] for x in range(len(CCREList)) if x % 2 == 0]
        CCRE_EndCoords = [CCREList[y] for y in range(len(CCREList)) if y % 2 != 0]
        KnownCCREs = [(x, y) for x, y in zip(CCRE_StartCoords, CCRE_EndCoords)]

        #Initialize __init__ and populate the variables with user defined vars
        Init = InitializeVars(
            List=[r"F:\R_dir\ATAC_midbrain_PND0\ENCFF893JNS.bam"], 
            Start=22618000, End=22634000, Strnd="+", Bool=True,
            Chromosome="chr2", BinSize=100, AverageFrames=False, Thresh=0.33,
            Export=r"F:\Test\Sample_", 
            CCREs=KnownCCREs, 
            AccountForOpenSpread=500
        )
        Sort = Sorters()
        Exprt = ExportClass()
        Stats = Statistics()
        
        ###################
        #Count selection - choose to process reads directly or by rounding them to the nearest bin
        ###################
        CountSelection = "direct"
        if CountSelection == "rounding":
            Process = CountReadsByGrouping().VarLoads(Init)
        elif CountSelection == "direct":
            Process = CountReadsDirectly().VarLoads(Init)
        ReadCountsCCRE = GetReadCountsPerCCRE()
        ReadCountsCCRE.VarLoads(Init, InputFiles = Process)
        
        ###################
        #Remove known CCRE regions
        ###################
        CCRE_Removed = RemoveKnownCCRE().VarLoads(Init, InputFiles = Process)
        
        ###################
        #Load Stats functions
        ###################
        Sts = StatisticalProcessing()
        Sts.StatisticalFunctions(Stats, Inputs = CCRE_Removed)
        
        ###################
        #Set Alignment
        ###################
        PerformAlign = "y"
        if PerformAlign == "y":
            Start = 22618300
            End = 22618500
            RefSeq = []
            from Align import aligner
            Align = aligner.AlignReads(
                AlignmentStart=Start, 
                AlignmentEnd=End, 
                ReferenceSequence=RefSeq,
                MainFrame=FrameList
                ) 
            AlignFxn = aligner.ReadAlignment()
            AlignFxn.ReadAlign(Align)
            
        elif PerformAlign == "n":
            pass
        else:
            raise(TypeError("Please enter either 'y' for yes or 'n' in the PerformAlign for alignments"))
        ###################
        #Reset working directory to original
        ###################
        if "CWD" in globals():
            ResetCWDToOld().ResetFxn(InputCWD=CWD)
    else:
        raise(NameError("Please input either 'useCLI' or 'useIDE' in the FileInput.py VARLOAD_METHOD variable"))


