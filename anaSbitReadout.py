#!/bin/env python
"""
    Script to analyze sbitreadout data
    By: Caterina Aruta (caterina.aruta@cern.ch) and Francesco Ivone (francesco.ivone@cern.ch)
    """
if __name__ == '__main__':
    import ROOT as r
    import numpy as np
    import root_numpy as rp
    import glob
    import sys
    import os
    import time
    from optparse import OptionParser
    from subprocess import call
    from array import array
    from gempython.utils.nesteddict import nesteddict as ndict
    from gempython.gemplotting.mapping.chamberInfo import chamber_iEta2VFATPos
    from gempython.gemplotting.mapping.chamberInfo import chamber_vfatPos2iEtaiPhi as vfat_to_etaphi
    from gempython.gemplotting.utils.anautilities import make3x8Canvas, saveSummaryByiEta, getMapping
    from gempython.gemplotting.utils.anaoptions import parser
    import pkg_resources

    parser.set_defaults(outfilename="sbitReadOut.root")
    (options, args) = parser.parse_args()
    path = options.filename
    size = ((options.GEBtype).lower())

    if not os.path.isdir(path):
        raise OSError(2, "No such file or directory", path)
    else:
        pass

    if size not in ('long', 'short'):
        raise AssertionError("Invalid value of GEBtype")
    else:
        pass

    filename = options.filename+"/sbitReadOut"
    outfilename = options.outfilename
    print("Analyzing: '%s'" % path)
    os.system("mkdir "+filename)

    """
    At the moment the output of the sbitReadout.py are .dat files with headers that will automatically fill the ROOT TTree. 
    However ther is one last ':' in the first line that shouldn't be there; consequently the ReadFile function is not able to understand the header. So before reading these .dat files, one has to be sure to remove the ':'
    To achieve this the following command is needed (removes the last : in the header of all .dat files in the input path)
    """
    os.system("find "+path+" -iname \*.dat -exec sed -i 's/7\/i:/7\/i/g' {} \;")

    # Set default histo behavior
    r.TH1.SetDefaultSumw2(False)
    r.gROOT.SetBatch(True)
    r.gStyle.SetOptStat(1111111)

    mappath = pkg_resources.resource_filename(
        'gempython.gemplotting', 'mapping/')
    # Loading the dictionary with the mapping
    vfat_ch_strips = getMapping(
        mappath+size+"ChannelMap_VFAT3-HV3b-V1_VFAT3-HV3b-V2.txt", False)
    print("\nVFAT channels to strips "+size+":\t MAP loaded")

    # Loading and reversing the dictionary with (eta , phi) <-> vfatN
    etaphi_to_vfat = ndict()
    for i in range(1, 9):
        etaphi_to_vfat[i] = {v: k for k,
                             v in chamber_iEta2VFATPos[i].iteritems()}

    """
        Now it's time to load all the input files and to merge them into one root TTree
        The TTree file it's going to auto extend each time a new file is found
        A TFile it's going to hold this TTree

        """
    # Creating the output File and TTree
    outF = r.TFile(filename+'/'+outfilename, 'recreate')
    inT = r.TTree('Packed', 'Tree Holding packed raw data')

    # searching for all the files with this format and adding them to the TTree
    start_time = time.time()
    print ("\nReading .dat files from the folder '%s'" % path)
    for idx, file in enumerate(glob.glob(path+'/sbitReadOut_run*.dat')):
        inT.ReadFile(file)
        # inT.Fill()

    print idx+1, 'input files have been read and added to the TTree'
    # Add the TTree into the TFile
    inT.Write()
    print 'TTree written\n'

    """
    Going to build the output tree starting from the previous TTree converted into an array.
    First of all, going to initilize the array which will hold the data

    """

    # copying the branch names in order to work with input TTree as an array
    branches = inT.GetListOfBranches()
    brName = np.empty([branches.GetEntries()], dtype='object')
    for i in range(0, branches.GetEntries()):
        brName[i] = branches[i].GetName()

    # converting the input tree in array then intialiting the unpackd TTree
    rawData = rp.tree2array(tree=outF.Packed, branches=brName)
    outT = r.TTree('unPacked', 'Tree holding unpacked data')

    """
    BRANCH VARIABLES DEFINITION
    vfatCH and strip have bigger size because they can hold different
    number of data: from 2 up to 16 (depending on cluster size). In this
    way the associated branch can be filled with these data in one line
    - 16 it's the max value allowed!
    - In each run it will be filled from 0 up to chHitPerCluster
    """
    vfatN = array('i', [0])
    chHitPerCluster = array('i', [0])
    vfatCH = array('i', 16*[0])
    strip = array('i', 16*[0])
    sbitSize = array('i', [0])
    L1Delay = array('i', [0])

    outT.Branch('vfatN', vfatN, 'vfatN/I')
    outT.Branch('chHitPerCluster', chHitPerCluster, 'chHitPerCluster/I')
    outT.Branch('vfatCH', vfatCH, 'vfatCH[chHitPerCluster]/I')
    outT.Branch('strip', strip, 'strip[chHitPerCluster]/I')
    outT.Branch('sbitSize', sbitSize, 'sbitSize/I')
    outT.Branch('L1Delay', L1Delay, 'L1Delay/I')

    """
    Defining both VFAT and iEta histos
    """
    # initializing vfat 1Dhisto
    # While strip & ch branch are filled with arrays, histos are filled with one entries at a time

    vfat_h_strip = ndict()
    vfat_h_ch = ndict()
    vfat_h_sbitSize = ndict()
    vfat_h_delay = ndict()

    for i in range(0, 24):
        vfat_h_strip[i] = r.TH1F("Strips_vs_hit{0}".format(
            i), "vfat{0}".format(i), 128, 0., 128.)
        vfat_h_ch[i] = r.TH1F("Chann_vs_hit{0}".format(
            i), "vfat{0}".format(i), 128, 0., 128.)
        vfat_h_sbitSize[i] = r.TH1F("Sbitsize_vs_hit{0}".format(
            i), "vfat{0} sbitSize".format(i), 8, 0., 8.)
        vfat_h_delay[i] = r.TH1F("L1A_Sbit_delay{0}".format(
            i), "vfat{0} L1A delay".format(i), 4096, 0., 4096.)
        vfat_h_ch[i].SetXTitle("Chann Num")
        vfat_h_strip[i].SetXTitle("Strip Num")
        vfat_h_strip[i].SetFillColorAlpha(r.kBlue, 0.35)
        vfat_h_ch[i].SetFillColorAlpha(r.kBlue, 0.35)

    # initializing eta 1Dhisto
    ieta_h_strip = ndict()
    ieta_h_ch = ndict()
    ieta_h_sbitSize = ndict()
    ieta_h_delay = ndict()

    for i in range(1, 9):
        ieta_h_strip[i] = r.TH1F("eta_Strips_vs_hit{0}".format(
            i), "i#eta = {0} | i#phi (1,2,3)".format(i), 384, 0., 384.)
        ieta_h_ch[i] = r.TH1F("eta_Chann_vs_hit{0}".format(
            i), "i#eta = {0} | i#phi (1,2,3)".format(i), 384, 0., 384.)
        ieta_h_sbitSize[i] = r.TH1F("eta_Sbitsize_vs_hit{0}".format(
            i), "i#eta = {0} sbitSize".format(i), 8, 0., 8.)
        ieta_h_delay[i] = r.TH1F("eta_L1A_Sbit_delay{0}".format(
            i), "i#eta = {0} L1A delay".format(i), 4096, 0., 4096.)
        ieta_h_strip[i].SetFillColorAlpha(r.kBlue, 0.35)
        ieta_h_ch[i].SetFillColorAlpha(r.kBlue, 0.35)
        ieta_h_strip[i].SetXTitle("Strip num")
        ieta_h_ch[i].SetXTitle("Chan num")

    # initializing 2Dhisto
    hh_ieta_strip = ndict()
    hh_ieta_ch = ndict()
    hh_ieta_strip[0] = r.TH2I(
        'ieta_strip', 'Strips summary        (i#phi = 1,2,3);strip number;i#eta', 384, 0, 384, 8, 0.5, 8.5)
    hh_ieta_ch[0] = r.TH2I(
        'ieta_ch', 'Channels summary        (i#phi = 1,2,3);chan number;i#eta', 384, 0, 384, 8, 0.5, 8.5)

    # loop over all branch names but the first (evnt num)
    for branch in brName[1:]:
        print "Unpacking now", branch + \
            ",", str(rawData[branch].shape).translate(None, "(),"), "entries"
        # h is a dummy variable to be printed in case of error
        for h, word in enumerate(rawData[branch]):

            sbitAddr = ((word) & 0x7FF)
            # INVALID ADDRESS CHECK
            if sbitAddr >= 1536:
                continue

            vfatN[0] = (7 - sbitAddr/192 + ((sbitAddr % 192)/64)*8)
            sbitSize[0] = ((word >> 11) & 0x7)
            chHitPerCluster[0] = 2*(sbitSize[0]+1)
            L1Delay[0] = ((word >> 14) & 0xFFF)
            eta = vfat_to_etaphi[vfatN[0]][0]
            phi = vfat_to_etaphi[vfatN[0]][1]
            # SBIT always includes doublet of adjacent channels
            vfatCH[0] = 2*(sbitAddr % 64)
            vfatCH[1] = vfatCH[0] + 1
            strip[0] = vfat_ch_strips[vfatN[0]]['Strip'][vfatCH[0]]
            strip[1] = vfat_ch_strips[vfatN[0]]['Strip'][vfatCH[1]]
            # In case of wrong mapping adjacent channels may not be adjacent strips, which is physically inconsistent
            if(np.abs(strip[0] - strip[1]) > 1):
                print 'WARNING: not adjacent strips'
                time.sleep(3)

            # filling vfat 1Dhistos
            vfat_h_strip[vfatN[0]].Fill(strip[0])
            vfat_h_strip[vfatN[0]].Fill(strip[1])
            vfat_h_ch[vfatN[0]].Fill(vfatCH[0])
            vfat_h_ch[vfatN[0]].Fill(vfatCH[1])
            vfat_h_delay[vfatN[0]].Fill(L1Delay[0])
            vfat_h_sbitSize[vfatN[0]].Fill(sbitSize[0])

            # filling ieta 1Dhistos
            ieta_h_strip[eta].Fill((phi-1)*128+strip[0])
            ieta_h_strip[eta].Fill((phi-1)*128+strip[1])
            ieta_h_ch[eta].Fill((phi-1)*128+vfatCH[0])
            ieta_h_ch[eta].Fill((phi-1)*128+vfatCH[1])
            ieta_h_delay[eta].Fill(L1Delay[0])
            ieta_h_sbitSize[eta].Fill(sbitSize[0])

            # filling 2Dhisto
            hh_ieta_strip[0].Fill(strip[0]+128*(phi-1), eta)
            hh_ieta_strip[0].Fill(strip[1]+128*(phi-1), eta)
            hh_ieta_ch[0].Fill(vfatCH[0]+128*(phi-1), eta)
            hh_ieta_ch[0].Fill(vfatCH[1]+128*(phi-1), eta)

            """
            A single sbit word can specify up to 16 adjacent channel hits
            The following loop takes care of this possibility
            """
            for i in range(2, chHitPerCluster[0]):
                # filling the adjacent channel
                vfatCH[i] = vfatCH[i-1] + 1
                # if the new channel exceeds the total VFAT channels, increase phi and move to the first cannel of the next VFAT
                if vfatCH[i] >= 128 and phi < 3:
                    phi = phi + 1
                    vfatCH[i] = 0
                    vfatN[0] = etaphi_to_vfat[eta][phi]
                # if the maximum of phi is reached (so there is no "next VFAT"), there must be some kind of error
                elif vfatCH[i] >= 128 and phi >= 3:
                    # ERROR
                    print "ERROR"
                    print "word", word
                    print "VFATN:", vfatN[0]
                    print "VFATCH:", vfatCH[0]
                    print "Cluster Sz:", sbitSize[0]
                    print "Eta: ", eta
                    print "phi: ", phi
                    break

                else:
                    pass
                # updating the strip
                strip[i] = vfat_ch_strips[vfatN[0]]['Strip'][vfatCH[i]]

                # At this point both strip and ch are updated, going to fill the histos
                vfat_h_strip[vfatN[0]].Fill(strip[i])
                vfat_h_ch[vfatN[0]].Fill(vfatCH[i])

                ieta_h_strip[eta].Fill((phi-1)*128+strip[i])
                ieta_h_ch[eta].Fill((phi-1)*128+vfatCH[i])
                hh_ieta_strip[0].Fill(strip[i]+128*(phi-1), eta)
                hh_ieta_ch[0].Fill(vfatCH[i]+128*(phi-1), eta)

            outT.Fill()

        pass

    #
    # Summaries Canvas
    #
    # make3x8Canvas
    canv_3x8 = make3x8Canvas(
        name="Strip_3x8canv",
        initialContent=vfat_h_strip,
        initialDrawOpt="hist",
        secondaryContent=None,
        secondaryDrawOpt="hist")
    canv_3x8.SaveAs(filename+'/StripSummary.png')

    canv_3x8 = make3x8Canvas(
        name="Chann_3x8canv",
        initialContent=vfat_h_ch,
        initialDrawOpt="hist",
        secondaryContent=None,
        secondaryDrawOpt="hist")
    canv_3x8.SaveAs(filename+'/ChannSummary.png')

    # Uncomment for more png plots
    """
    canv_3x8 = make3x8Canvas(
                             name="SbitSize_3x8canv",
                             initialContent=vfat_h_sbitSize,
                             initialDrawOpt="hist",
                             secondaryContent=None,
                             secondaryDrawOpt="hist")
    canv_3x8.SaveAs(filename+'/SbitSizeSummary.png')

    canv_3x8 = make3x8Canvas(
                             name="L1A_Delay_3x8canv",
                             initialContent=vfat_h_delay,
                             initialDrawOpt="hist",
                             secondaryContent=None,
                             secondaryDrawOpt="hist")
    canv_3x8.SaveAs(filename+'/L1A_DelaySummary.png')
    """

    # saveSummaryByiEta
    saveSummaryByiEta(ieta_h_strip, name='%s/ietaStripSummary.png' %
                      filename, trimPt=None, drawOpt="")
    saveSummaryByiEta(ieta_h_ch, name='%s/ietaChanSummary.png' %
                      filename, trimPt=None, drawOpt="")
    saveSummaryByiEta(ieta_h_sbitSize, name='%s/ietaSbitSizeSummary.png' %
                      filename, trimPt=None, drawOpt="")
    saveSummaryByiEta(ieta_h_delay, name='%s/ietaDelaySummary.png' %
                      filename, trimPt=None, drawOpt="")

    # Making&Filling folders in the TFile
    outT.Write()
    vfatDir = outF.mkdir("VFATs")
    ietaDir = outF.mkdir("iETAs")

    vfatDir.cd()
    for vfat in range(0, 24):
        tempDir = vfatDir.mkdir("VFAT%i" % vfat)
        tempDir.cd()
        vfat_h_strip[vfat].Write()
        vfat_h_ch[vfat].Write()
        vfat_h_delay[vfat].Write()
        vfat_h_sbitSize[vfat].Write()

    ietaDir.cd()
    for ieta in range(1, 9):
        tempDir = ietaDir.mkdir("iETA%i" % ieta)
        tempDir.cd()
        ieta_h_strip[ieta].Write()
        ieta_h_ch[ieta].Write()
        ieta_h_delay[ieta].Write()
        ieta_h_sbitSize[ieta].Write()

    # 2D histos aesthetics
    line1 = r.TLine(128, 0.5, 128, 8.5)
    line2 = r.TLine(256, 0.5, 256, 8.5)
    line1.SetLineColor(r.kRed)
    line1.SetLineWidth(3)
    line2.SetLineColor(r.kRed)
    line2.SetLineWidth(3)

    canv = r.TCanvas("summary", "summary", 500*8, 500*3)
    canv.SetGridy()
    canv.cd()
    hh_ieta_strip[0].Draw('9COLZ')
    line1.Draw()
    line2.Draw()
    canv.Update()
    canv.SaveAs(filename+'/StripvsiEta.png')

    canv.Clear()
    canv.cd()
    hh_ieta_ch[0].Draw('COLZ')
    line1.Draw()
    line2.Draw()
    canv.Update()
    canv.SaveAs(filename+'/ChvsiEta.png')

    outF.Close()
    print "\n---Took", (time.time() - start_time) / \
        int(idx), "seconds for each .dat file---"
    print "\nGaranting permission to ", filename, "..."
    os.system("chmod -R 770 "+filename)
    print "Data stored in", filename+'/'+outfilename, "\nBye now"
