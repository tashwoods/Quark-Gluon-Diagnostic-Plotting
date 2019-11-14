#!/usr/bin/python
import ROOT, os, math, glob, sys
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
from ROOT import TCanvas, TColor, TGaxis, TH1F, TPad, TFile, TPaveText, TPaveLabel
#import numpy as np
from ROOT import THStack
from ROOT import gROOT
from ROOT import TPad
#import matplotlib.pyplot as plt
from array import array
from filefunctions import *
from classes import *
import tempfile
import numpy as np
from decimal import Decimal

def yieldfunctions():
    #Specify slope and intercept used to distinguish quarks and gluons via nTrk and signal file to run on
    slope = sys.argv[1]
    inter = sys.argv[2]
    signal = ["1lep_HVTWW"]

    #Set up configuration for output files
    outputdir = "output"
    check_if_directory_exists(outputdir)
    eventfile = outputdir + '/eventyields.txt'
    eventtxtfile = open(eventfile,"w")
    eventtxtfile.write('\t\t\t\t\t\t\t\t\t\t\t'+'QG Tag Event Yield\n')
    jetfile = outputdir + '/jetyields.txt'
    jettxtfile = open(jetfile,"w")
    jettxtfile.write('\t\t\t\t\t\t'+'%taggable jets\t%Q jets\t\t%Q jets Q-tag\t\t%Q jets G-tag\t\t%G jets\t\t\t%G jets G-tag\t\t%G jets Q-tag\n')
    hfilepath = outputdir + "/histos.root"
    hfile = TFile( hfilepath, 'RECREATE', 'Root file with my histograms')

    #Do the real work! (making plots)
    for s in signal:
        sig = signalObject(s, slope, inter)
        make_standard_plots(sig,hfile,outputdir, eventtxtfile, jettxtfile)
        del sig

    #cleanup
    hfile.Write()
    hfile.Close()
    jettxtfile.close()
    eventtxtfile.close()
    del hfile

def make_standard_plots(sig, hfile, outputdir, eventtxtfile, jettxtfile):
    #make 1d variable plots
    for v in sig.variables:
        #makes simple overlay of input files for given variable
        diagnostic_1dplots(v, sig.files, hfile, outputdir, sig) 
        #plots overlays with QG requirements and ratio of QG/Nominal
        sophisticated_1dplots(v, sig.files, hfile, outputdir, sig) 

    #loop over signal and background files
    for f,filename in zip(sig.files, sig.shortfilenames):
        #obtain jet yields and nTrk vs pT plots per jet
        jettxtfile.write(filename + '\n')
        for j in range(len(sig.ntrkvar)):
            #makes jet text file with jet composition and tag/mistag efficiencies 
            make_jet_table(f, sig, j, hfile, outputdir, jettxtfile) 
            #makes 2-d variables plots (nominal: pT vs nTrk) for quarks, gluons, and all particles
            make_ntrk_var_plots(f,j,sig,sig.heatmapvar, hfile, outputdir, "quark") 
            make_ntrk_var_plots(f,j,sig,sig.heatmapvar, hfile, outputdir, "gluon") 
            make_ntrk_var_plots(f,j,sig,sig.heatmapvar, hfile, outputdir, "all") 

        #make event yield table
        eventtxtfile.write("%20s" %filename)
        make_event_table(f, sig, sig.eventvar, hfile, outputdir, eventtxtfile)

    #Plot m(lvqq) distribution with tagging, and its impact on S/sqrt(b)
    make_evtvar_plots(hfile, outputdir, sig)

def make_evtvar_plots(hfile, outputdir, sig):
    #create canvas, legend, and temp file to store them and histograms created later
    tempf = ROOT.TFile.Open(tempfile.mktemp(), 'RECREATE')
    c2, pad1, pad2, pad3 = create_canvas_three_pads(tempf)
    pad1.cd()
    mleg = ROOT.TLegend(0.35,0.8,0.9,0.9)
    mleg.SetBorderSize(0)

    #create matrix to store histograms created in loop
    np_hist_array = np.empty([len(sig.eventqgselection), len(sig.files)], dtype=object)
    

    #iterate over QG Selections
    for j in range(len(sig.eventqgselection)):
        #iterate over signal and background files
        for myfile,filename,color,thisindex in zip(sig.files, sig.shortfilenames, sig.color, sig.indices):
            #obtain tree for given file
            f = ROOT.TFile.Open(myfile, 'r')
            t = f.Get(sig.treename)
            
            if j ==0: #Nominal Selection
                selection=sig.weight + "*(" + sig.selection + ")"
            else: #QG Selection(s)
                selection=sig.weight + "*(" + sig.selection + ")*(" + sig.eventqgselection[j] + ")"

            histname=filename #to keep track of histname without variable name duplicates
            if 'tag' in sig.namesqgselection[j]:
                histname += " QG Truth Tag"
            hist1 = ROOT.TH1F(histname,"", int(len(sig.binning)-1), array('d',sig.binning))
            
            #Plot Aesthetics
            if j==1: #Nominal Plot uses solid line
                hist1.SetLineStyle(2)
            if j==2: #QG Selection use dashed lines
                hist1.SetLineStyle(9)
            hist1.SetLineColor(color)
            hist1.SetLineWidth(2)
            hist1.SetLineColorAlpha(color, 1-(0.3*j))
            hist1.SetTitle(sig.eventqgselection[j])
            hist1.SetStats(0)
            hist1.GetXaxis().SetTitle(get_var_title(sig.eventvar))
            hist1.GetYaxis().SetTitleOffset(0.1)
            hist1.GetYaxis().SetTitleSize(0.1)

            #Use tree inputs to make histograms
            makehisto = sig.eventvar + ">>" + histname
            hist1.Scale(sig.filesxs[j])
            t.Draw(makehisto, selection, "goff")
            hist1.SetDirectory(tempf)

            #save histograms to array for S/sqrt(B) and QGTag/Nominal calculations below
            np_hist_array[j][thisindex] = hist1
            f.Close()
            del hist1, f, t, histname, makehisto

    #create stacked distribution of signals + backgrounds and draw on pad 1
    stack_hist_array, mleg = create_stacked_histo_array(np_hist_array)
    pad1.SetLogy()
    draw_histo_from_matrix(stack_hist_array, "Number of Weighted Events")
    mleg.Draw()

    #calculate QGTag/Nominal histograms and plot on pad 2
    pad2.cd()
    ratio_array = ratio_hist_from_matrix(np_hist_array)
    draw_histo_from_matrix(ratio_array, "QG Tag/Nominal")

    #calculate the change in S/sqrt(B) with QG Tagging and plot on pad 3
    pad3.cd()
    s_root_b_matrix = get_s_root_b_matrix(np_hist_array)
    ratio_s_root_b = ratio_hist_from_matrix(s_root_b_matrix)
    draw_histo_from_matrix(ratio_s_root_b, "\Delta (S/\sqrt{B})")

    #Save Plot
    c2.Print(outputdir + "/" + sig.signal + "_" + sig.eventvar + "event.pdf")
    hfile.cd()
    c2.Write(sig.signal + "_" + sig.eventvar)

    return 

def get_s_root_b_matrix(mymatrix):
    np_hist_array = np.empty([mymatrix.shape[0], mymatrix.shape[1]-1], dtype=object)
    for i in range(mymatrix.shape[0]):
        #calculate sqrt(B) histogram
        sqrtbkg = sqrthisto(mymatrix[i][mymatrix.shape[1]-1])
        #iterate over signals (not background as well hence the "-1")
        for j in range(mymatrix.shape[1]-1): 
            #divide signal histogram (mymatrix[i][j]) by sqrt(bkg) histgram
            thishist = get_ratio_histo(sqrtbkg, mymatrix[i][j], mymatrix[i][j].GetName()) 
            np_hist_array[i][j]=thishist
            del thishist
    return np_hist_array

def draw_histo_from_matrix(mymatrix, yaxistitle):
    thismax = np_get_max(mymatrix)
    for i in range(mymatrix.shape[0]):
        for j in range(mymatrix.shape[1]):
            if i ==0 and j ==0:
                mymatrix[0][0].SetMaximum(1.1*thismax)
                if thismax < 2:
                    mymatrix[0][0].SetMinimum(0.5)
                mymatrix[i][j].Draw("HIST")
                mymatrix[i][j].GetXaxis().SetLabelSize(0.1)
                mymatrix[i][j].GetYaxis().SetLabelSize(0.08)
                mymatrix[i][j].GetXaxis().SetTitleOffset(1.2)
                mymatrix[i][j].GetXaxis().SetTitleSize(0.1)
                mymatrix[i][j].GetYaxis().SetTitle(yaxistitle)
                mymatrix[i][j].GetYaxis().SetTitleOffset(0.3)
                mymatrix[i][j].GetYaxis().SetTitleSize(0.1)
            else:
                mymatrix[i][j].Draw("HIST SAME")
    return

def draw_markers_from_matrix(mymatrix):
    for i in range(mymatrix.shape[0]):
        for j in range(mymatrix.shape[1]):
            if i ==0 and j ==0:
                mymatrix[i][j].Draw()
            else:
                mymatrix[i][j].Draw("SAME")


def create_stacked_histo_array(mymatrix):
    #creating 1d array of stacked histograms for each analysis setup (e.g. number of rows in mymatrix)
    stack_hist_array = np.empty([mymatrix.shape[0], mymatrix.shape[1]-1], dtype=object) 

    mleg = ROOT.TLegend(0.65,0.8,0.75,0.99)
    mleg.SetBorderSize(0)
    for i in range(mymatrix.shape[0]): #iterate over variations (rows)
        for j in range(mymatrix.shape[1]-1): #iterate over signal files to draw
            hs = THStack(str(mymatrix[i][j].GetName()), "") #create stack histo for each variation
            hs.Add(mymatrix[i][mymatrix.shape[1]-1]) #add background first
            hs.Add(mymatrix[i][j])
            if j ==0:
                #only draw background once since it is independent of signal sample
                mleg.AddEntry(mymatrix[i][mymatrix.shape[1]-1], mymatrix[i][mymatrix.shape[1]-1].GetName(), "L") 
            mleg.AddEntry(mymatrix[i][j],mymatrix[i][j].GetName(), "L")
            stack_hist_array[i][j]=hs
            del hs
    return stack_hist_array, mleg

def ratio_hist_from_matrix(mymatrix):
    #create empty histogram that will be filled with ratio of input histograms
    ratio_array = np.empty([mymatrix.shape[0], mymatrix.shape[1]], dtype=object)
    for i in range(mymatrix.shape[1]):
        for j in range(mymatrix.shape[0]):
            #calculate ratio of nominal histogram (mymatrix[0][i] to mymatrix[j][i] and save in ratio_array
            ratio_array[j][i]= get_ratio_histo(mymatrix[0][i], mymatrix[j][i], mymatrix[j][i].GetName())
    return ratio_array
    
def printhistobins(hist):
    for i in range(hist.GetNbinsX()):
        print 'bin ' + str(i) + ": " + str(hist.GetBinContent(i))

def getratio(hist1, hist2):
    #returns histogram that is hist2/hist1
    rhist=hist1.Clone("rhist")
    rhist.SetStats(0)
    rhist.Divide(hist2)
    return rhist

def sqrthisto(hist):
    #returns histogram that is the square root of input histogram
    hist2=hist.Clone("hist2")
    for i in range(hist2.GetNbinsX()+1):
        if hist2.GetBinContent(i) > 0:
            hist2.SetBinContent(i,math.sqrt(hist2.GetBinContent(i)))
        else:
            hist2.SetBinContent(i,0)
    hist2.SetBinContent(hist2.GetNbinsX()+2, 0)
    return hist2

def make_event_table(myfile, sig, var, hfile, outputdir, txtfile):
    #define selections
    preselection = sig.weight + "*(" + sig.selection + ")*"
    qgselections = [preselection + s for s in sig.eventqgselection]
    for x in range(len(qgselections)):
        if qgselections[x][-1] == '*':
            qgselections[x] = qgselections[x][:-1]

    #obtain trees from file
    index1 = myfile.rfind("/") + 1
    index2 = myfile.rfind(".")
    filename = myfile[index1:index2]
    f = ROOT.TFile.Open(myfile, 'r')
    t = f.Get(sig.treename)

    tempf = ROOT.TFile.Open("temp.root", 'RECREATE')
    mybins, mymin, mymax = get_histo_params(var)
    weights_dict = {}

    selections_dict = {
        "nominal": qgselections[0],
        "QG Tag Sig Jets": qgselections[1]
    }

    if 'VBF' in sig.signal:
        selections_dict["QG Tag Sig and VBF jets"] = qgselections[2]

    for x in selections_dict:
        hist = ROOT.TH1F("hist","", mybins, mymin, mymax)
        makehisto = var + ">>hist"
        t.Draw(makehisto, selections_dict[x], "goff")
        hist.SetDirectory(tempf)
        weights_dict[x]=hist.GetSumOfWeights()
        del hist

    #Tagging Signal Jets yields:
    yield_tag_sig_jets = 100*(weights_dict["QG Tag Sig Jets"]/weights_dict["nominal"])

    #Tagging Signal and VBF Jet Yields (only for VBF signals)
    if 'VBF' in sig.signal:
        yield_tag_sig_vbf_jets = 100*(weights_dict["QG Tag Sig and VBF jets"]/weights_dict["nominal"])
        txtfile.write('\t\t\t%.2E\t\t\t%.2E\t\t\t%4.2f\t\t\t%.2E\t\t\t%4.2f\n' % (weights_dict["nominal"],weights_dict["QG Tag Sig Jets"], yield_tag_sig_jets, weights_dict["QG Tag Sig and VBF jets"],yield_tag_sig_vbf_jets))
    else:
        txtfile.write('\t\t\t\t\t%4.2f\n' % (yield_tag_sig_jets))
    return 

def make_jet_table(myfile, sig, j, hfile, outputdir, jettxtfile):
    #define selections
    preselection = sig.weight + "*(" + sig.selection + ")"
    if j==0:
        taggable_selection = preselection + "*(!" + sig.nottaggable1 + ")"
    if j==1:
        taggable_selection = preselection + "*(!" + sig.nottaggable2 + ")"

    #truth definitions
    quark_selection = "*(" + sig.truthvar[j] + " > 0  && " + sig.truthvar[j] + " < 5 )"
    gluon_selection = "*(" + sig.truthvar[j] + "==21)"
    taggable_quark_selection = taggable_selection + quark_selection
    taggable_gluon_selection = taggable_selection + gluon_selection

    if j==0:
        quarktag = "*(" + sig.sigJ1_rqmt + ")"
        gluontag = "*(!" + sig.sigJ1_rqmt + ")"
    if j==1:
        quarktag = "*(" + sig.sigJ2_rqmt + ")"
        gluontag = "*(!" + sig.sigJ2_rqmt + ")"


    taggable_quark_quark_tagged = taggable_quark_selection + quarktag
    taggable_quark_gluon_tagged = taggable_quark_selection + gluontag
    taggable_gluon_quark_tagged = taggable_gluon_selection + quarktag
    taggable_gluon_gluon_tagged = taggable_gluon_selection + gluontag

    selections_dict = {
        "SR_selection": preselection,
        "taggable": taggable_selection,
        "taggable quarks": taggable_quark_selection,
        "taggable gluons": taggable_gluon_selection,
        "Quarks tagged quarks": taggable_quark_quark_tagged,
        "Quarks tagged gluons": taggable_quark_gluon_tagged,
        "Gluons tagged quarks": taggable_gluon_quark_tagged,
        "Gluons tagged gluons": taggable_gluon_gluon_tagged
    }

    #obtain trees from file
    index1 = myfile.rfind("/") + 1
    index2 = myfile.rfind(".")
    filename = myfile[index1:index2]
    f = ROOT.TFile.Open(myfile, 'r')
    t = f.Get(sig.treename)

    tempf = ROOT.TFile.Open("temp.root", 'RECREATE')
    mybins, mymin, mymax = get_histo_params(sig.eventvar)
    weights_dict = {}

    for x in selections_dict:
        hist = ROOT.TH1F("hist","", mybins, mymin, mymax)
        makehisto = sig.eventvar + ">>hist"
        t.Draw(makehisto, selections_dict[x], "goff")
        hist.SetDirectory(tempf)
        weights_dict[x]=hist.GetSumOfWeights()
        del hist

    #Taggable Jets
    percent_jets_taggable = 100*(weights_dict["taggable"]/weights_dict["SR_selection"])

    #True quark jets
    percent_jets_quarks_taggable = 100*(weights_dict["taggable quarks"]/weights_dict["taggable"])
    percent_quarks_jets_quark_tagged = 100*(weights_dict["Quarks tagged quarks"]/weights_dict["taggable quarks"])
    percent_quarks_jets_gluon_tagged = 100*(weights_dict["Quarks tagged gluons"]/weights_dict["taggable quarks"])

    #True gluon jets
    percent_jets_gluons_taggable = 100*(weights_dict["taggable gluons"]/weights_dict["taggable"])
    if weights_dict["taggable gluons"]!=0:
        percent_gluon_jets_gluon_tagged = 100*(weights_dict["Gluons tagged gluons"]/weights_dict["taggable gluons"])
        percent_gluon_jets_quark_tagged = 100*(weights_dict["Gluons tagged quarks"]/weights_dict["taggable gluons"])
    else:
        percent_gluon_jets_gluon_tagged = 0
        percent_gluon_jets_quark_tagged = 0
    jettxtfile.write(sig.jets[j]+'\t\t\t\t\t%4.2f\t\t\t\t%4.2f\t\t\t%4.2f\t\t\t\t\t\t%4.2f\t\t\t\t\t\t%4.2f\t\t\t\t\t%4.2f\t\t\t\t\t%4.2f\n' % (percent_jets_taggable, percent_jets_quarks_taggable, percent_quarks_jets_quark_tagged, percent_quarks_jets_gluon_tagged, percent_jets_gluons_taggable, percent_gluon_jets_gluon_tagged, percent_gluon_jets_quark_tagged))
    return 

def make_ntrk_var_plots(myfile, j, sig, var, hfile, outputdir, particles):
    fit2dhisto=True
    #define selections
    if particles=='quark':
      selection = sig.weight + "*(" + sig.truthvar[j] + ">0)*(" + sig.truthvar[j] + "<6)*(" + sig.selection + ")"
    elif particles=='gluon':
      selection = sig.weight + "*(" + sig.truthvar[j] + "==21)*(" + sig.selection + ")"
    elif particles=='all':
      selection = sig.weight + "*(" + sig.selection + ")"
    selectionpass = selection + "*(" + sig.ntrkvar[j] + " > " + sig.ntrkmin + " && " + sig.ntrkvar[j] + "<= (" + str(sig.qgslope) + "*TMath::Log(" + sig.ptvar[j] + ") + " + str(sig.qginter)+ "))"

    #obtain trees from file
    index1 = myfile.rfind("/") + 1
    index2 = myfile.rfind(".")
    filename = myfile[index1:index2]
    f = ROOT.TFile.Open(myfile, 'r')
    t = f.Get(sig.treename)

    tempf = ROOT.TFile.Open("temp.root", 'RECREATE')

    qgstring = str(sig.qgslope) + '*x + ' + str(sig.qginter) 
    f1 = ROOT.TF1("f1",qgstring,0.,10.);
    f1.SetLineColor(4)

    #create nTrk pt TH2F
    c3 = ROOT.TCanvas("c3", "c3", 0,  0, 1500, 1000)
    logptbins, logptmin, logptmax = get_histo_params('logpt')
    ntrkbins, ntrkmin, ntrkmax = get_histo_params('nTrk')
    if var=="pt" or var=="e":
        variables = sig.ntrkvar[j] + ":TMath::Log(" + sig.ptvar[j] +  ")"
    if var == "eta":
        variables = 'FIX ME'

    hist1 = ROOT.TH2F("hist1","", logptbins, logptmin, logptmax, ntrkbins, ntrkmin, ntrkmax)
    makehisto = variables + ">>hist1"
    t.Draw(makehisto, selection, "goff")
    hist1.SetDirectory(tempf)

    #create nTrk pt TH2F for jets that pass cut
    histpass = ROOT.TH2F("histpass","", logptbins, logptmin, logptmax, ntrkbins, ntrkmin, ntrkmax)
    makehisto = variables + ">>histpass"
    t.Draw(makehisto, selectionpass, "goff")
    histpass.SetDirectory(tempf)
    del f,t

    #find QG cut yield
    npass = histpass.GetSumOfWeights()
    nall = hist1.GetSumOfWeights()
    if nall == 0:
        myyield = 0
    else:
        myyield = 100*(npass/nall)
    syield = '%.2f'%myyield
    legendname = 'QG Cut, yield = ' + str(syield) + "%"

    #Find Best fit line and store in f2 to it is compatiable with JSRoot
    if fit2dhisto:
        myfit = hist1.Fit("pol1", "S")
        inter = myfit.Parameter(0)
        slope = myfit.Parameter(1)
        tinter = '%.1f'%inter
        tslope = '%.1f'%slope
        mystring = "Best Fit (nTrk = " + str(tslope) + "ln(p_{T}) + " + str(tinter)  + ")"
        f2 = ROOT.TF1("f2", "[0] + [1]*x + [2]*x + [3]*x", 0, 10)
        f2.SetParameter(0, myfit.Parameter(0))
        f2.SetParameter(1, myfit.Parameter(1))


    #plotting aesthetics!
    canvasname = filename + "_" + sig.selection + "_" + sig.ntrkvar[j]
    hist1.SetTitle(get_legend_hist_name(canvasname))
    hist1.Draw("COLZ")
    if fit2dhisto:
        f2.SetLineStyle(9)
        f2.SetLineColor(6)
        f2.Draw("same")
    leg = ROOT.TLegend(0.1,0.75,0.3,0.9)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.04)
    leg.AddEntry(f1,legendname, "l")
    if fit2dhisto:
        leg.AddEntry(f2,mystring, "l")
    leg.Draw()
    f1.Draw("same")
    if var == "pt":
        axisname = "ln(p_{T}) [GeV]"
    if var == "eta":
        axisname ="eta"
    hist1.GetXaxis().SetTitle(axisname)
    if 'J1' in sig.ntrkvar[j]:
        hist1.GetYaxis().SetTitle("Leading Jet Number of Reconstructed Tracks")
    if 'J2' in sig.ntrkvar[j]:
        hist1.GetYaxis().SetTitle("Sub-Leading Jet Number of Reconstructed Tracks")
    hist1.SetStats(0)
    hist1.Scale(1/hist1.GetSumOfWeights())
    c3.SetRightMargin(0.2)
    c3.Print(outputdir + "/" + canvasname + particles + ".pdf")
    hfile.cd()
    c3.Write(canvasname)
    return hist1

def diagnostic_1dplots(myvar, files, hfile, outputdir, sig):
    inselection = sig.selection
    treename = sig.treename
    weightname = sig.weight
    qgselection = sig.eventqgselection

    #create canvas and legend 
    tempf = ROOT.TFile.Open(tempfile.mktemp(), 'RECREATE')
    c1 = ROOT.TCanvas("c1", "c1", 0,  0, 1000, 700)
    c1.SetLogy()

    mybins, mymin, mymax = get_histo_params(myvar)
    histo_array=[]
    yields=[]

    #iterate over signal and background files
    for i in range(len(files)):
        myfile = files[i]
        index1 = myfile.rfind("/") + 1
        index2 = myfile.rfind(".")
        filename = myfile[index1:index2]
        f = ROOT.TFile.Open(myfile, 'r')
        t = f.Get(treename)
        selection=weightname + "*(" + inselection + ")"
        histname = filename

        if 'X_resolved_WW_m' in myvar:
            hist1 = ROOT.TH1F(histname,"", int(len(sig.binning)-1), array('d',sig.binning))
        else:
            hist1 = ROOT.TH1F(histname,"", mybins, mymin, mymax)
            makehisto = myvar + ">>" + histname
            hist1.GetYaxis().SetTitle("Number of Weighted Events")
            t.Draw(makehisto, selection, "goff")
            
            hist1.SetLineColor(sig.color[i])
            hist1.SetLineWidth(2)
            hist1.SetDirectory(tempf)
            hist1.SetStats(0)
            hist1.GetXaxis().SetTitle(get_var_title(myvar))
            histo_array.append(hist1)
        f.Close()

    legxmin, legymin, legxmax, legymax = get_legend_values(myvar)
    mleg = ROOT.TLegend(legxmin, legymin, legxmax, legymax)
    mleg.SetBorderSize(0)
    
    #overlay plots
    c1.cd()
    mymax = get_max(histo_array)*1.3
    histo_array[i].StatOverflows(1)
    for i in range(len(histo_array)):
        if i==0:
            histo_array[i].SetMaximum(mymax)
            histo_array[i].Draw("HIST")
        else:
            histo_array[i].Draw("HIST SAME")
        histmean = histo_array[i].GetMean()
        histerror = histo_array[i].GetStdDev()
        thishistname = get_legend_hist_name(histo_array[i].GetName())
        if 'pdgid' not in myvar:
            thishistname +=' \mu = %5.2f +/- %5.2f' % (histmean, histerror)
        mleg.AddEntry(histo_array[i],thishistname, "L")

    mleg.Draw()
    c1.Update()

    c1.Print(outputdir +"/" + myvar+"_"+inselection+"_diag.pdf")
    c1.Write(myvar+"_"+inselection)

    hfile.cd()

    tempf.Write()
    tempf.Close()
    del c1, hist1, histo_array, mleg, tempf, f
    return



def sophisticated_1dplots(myvar, files, hfile, outputdir, sig):
    #rename object attributes locally to simplify
    inselection = sig.selection
    treename = sig.treename
    weightname = sig.weight
    qgselection = sig.eventqgselection

    #create temp file for objects, canvas, legend, object arrays
    tempf = ROOT.TFile.Open(tempfile.mktemp(), 'RECREATE')
    c2, pad1, pad2 = create_canvas_two_pads(tempf)
    pad1.cd()
    mybins, mymin, mymax = get_histo_params(myvar)
    histo_array=[]
    ratio_histo_array = []
    yields=[]

    #iterate over signal and background files
    for i in range(len(files)):
        myfile = files[i]
        index1 = myfile.rfind("/") + 1
        index2 = myfile.rfind(".")
        filename = myfile[index1:index2]
        f = ROOT.TFile.Open(myfile, 'r')
        t = f.Get(treename)
        for j in range(len(qgselection)): #iterate over qg selections
            if j ==0: #first histogram is the nominal, so no QG requirement
                selection=weightname + "*(" + inselection + ")"
            else:
                selection=weightname + "*(" + inselection + ")*(" + qgselection[j] + ")"

            histname = filename + "_" + sig.namesqgselection[j]

            #determine histogram binning based on variable
            if 'X_resolved_WW_m' in myvar:
                hist1 = ROOT.TH1F(histname,"", int(len(sig.binning)-1), array('d',sig.binning))
            else:
                hist1 = ROOT.TH1F(histname,"", mybins, mymin, mymax)
            
            makehisto = myvar + ">>" + histname
            t.Draw(makehisto, selection, "goff")

            #Plotting Aesthetics
            hist1.SetLineColor(sig.color[i])
            hist1.SetLineWidth(2)
            if j ==1: #nominal distribution solid
                hist1.SetLineStyle(2)
            if j==2: #QG distribution dashed
                hist1.SetLineStyle(3)

            hist1.SetDirectory(tempf)
            hist1.SetTitle(inselection)
            hist1.SetStats(0)
            hist1.GetXaxis().SetTitle(myvar)
            hist1.GetYaxis().SetTitle("Number Of Weighted Events")
            histo_array.append(hist1)

            #Obtain ratio hist for QG selections, which will be plotted as a separate panel
            if j != 0:
                #divide my first histogram in array that is nominal-->i*len(qgselection) because i iterates over files
                ratio_hist = get_ratio_histo(histo_array[i*len(qgselection)], hist1, hist1.GetName()) 
                ratio_hist.SetDirectory(tempf)
                ratio_hist.GetXaxis().SetTitle(myvar)
                ratio_hist.GetYaxis().SetTitle("QG Tag/Nominal")
                ratio_hist.GetXaxis().SetLabelSize(0.1)
                ratio_hist.GetYaxis().SetLabelSize(0.1)
                ratio_hist.GetXaxis().SetTitleOffset(1.7)
                ratio_hist.GetXaxis().SetTitleSize(0.1)
                ratio_hist.GetYaxis().SetTitleOffset(0.32)
                ratio_hist.GetYaxis().SetTitleSize(0.1)
                ratio_hist.GetYaxis().SetRangeUser(0.5, 1.1)
                ratio_hist.SetLineColor(sig.color[i])
                ratio_histo_array.append(ratio_hist)
                if hist1.GetSumOfWeights() != 0:
                    yields.append(histo_array[i*len(qgselection)].GetSumOfWeights()/hist1.GetSumOfWeights())
                else:
                    yields.append(0)
        
        f.Close()

    #Make a legend and plot to keep track of distributions
    legxmin, legymin, legxmax, legymax = get_legend_values(myvar)
    mleg = ROOT.TLegend(legxmin, legymin, legxmax, legymax)
    mleg.SetBorderSize(0)
    
    #overlay nominal plots
    pad1.cd()
    mymax = get_max(histo_array)*1.1
    histo_array[i].StatOverflows(1)
    for i in range(len(histo_array)):
        if i==0:
            histo_array[i].SetMaximum(mymax)
            histo_array[i].Draw("HIST")
        else:
            histo_array[i].Draw("HIST SAME")
        histmean = histo_array[i].GetMean()
        histerror = histo_array[i].GetStdDev()
        thishistname = get_legend_hist_name(histo_array[i].GetName())
        thishistname +=' mean:%5.4f +/- %5.4f' % (histmean, histerror)
        mleg.AddEntry(histo_array[i],thishistname, "L")

    mleg.Draw()
    pad1.SetLogy()
    c2.Update()
    pad2.cd()

    #overlay ratio plots
    ratio_mymax = get_max(ratio_histo_array) + 0.01
    for myx in range(len(ratio_histo_array)):
        if myx==0:
            ratio_histo_array[myx].GetYaxis().SetRangeUser(-0.1, 1.1)
            ratio_histo_array[myx].Draw("")
        else:
            ratio_histo_array[myx].Draw("SAME")
    c2.Print(outputdir +"/" + myvar+"_"+inselection+".pdf")
    c2.Write(myvar+"_"+inselection)
    hfile.cd()

    tempf.Write()
    tempf.Close()
    del c2, pad1, pad2, hist1, ratio_hist, histo_array, ratio_histo_array, mleg, tempf, f
    return

def get_ratio_histo(hist1, hist2, histname): #this divides hist2/hist1
    numhist=hist2.Clone(histname)
    numhist.Divide(hist1)
    numhist.SetStats(0)
    numhist.SetTitle("")
    for i in range(numhist.GetNbinsX()+1):
        numhist.SetBinError(i, 0.1)
    return numhist

def get_single_yield_percentage(cut_sum, nom_sum):
    if nom_sum !=0: 
        myyield = 100*(cut_sum/nom_sum)
    else:
        myyield = 0
    return myyield

def create_canvas_two_pads(fp):
    fp.cd()
    c1 = ROOT.TCanvas("c1", "c1", 0,  0, 1000, 700)
    pad1 = TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0)
    pad1.SetGridx()
    pad1.Draw()
    c1.cd()  # returns to main canvas before defining pad2
    pad2 = TPad("pad2", "pad2", 0, 0.01, 1, 0.3)
    pad2.SetTopMargin(0)  # joins upper and lower plot
    pad2.SetBottomMargin(0.35)
    pad2.SetGridx()
    pad2.Draw()
    pad1.cd()
    return c1, pad1, pad2

def no_temp_file_create_canvas_two_pads():
    c1 = ROOT.TCanvas("c1", "c1", 0,  0, 1000, 700)
    pad1 = TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0)
    pad1.SetGridx()
    pad1.Draw()
    c1.cd()  # returns to main canvas before defining pad2
    pad2 = TPad("pad2", "pad2", 0, 0.01, 1, 0.3)
    pad2.SetTopMargin(0)  # joins upper and lower plot
    pad2.SetBottomMargin(0.35)
    pad2.SetGridx()
    pad2.Draw()
    pad1.cd()
    return c1, pad1, pad2

def create_canvas_three_pads(fp):
    fp.cd()
    c1 = ROOT.TCanvas("c1", "c1", 0,  0, 1000, 700)
    pad1 = TPad("pad1", "pad1", 0, 0.7, 1, 1.0)
    pad1.SetBottomMargin(0)
    pad1.SetGridx()
    pad1.Draw()
    c1.cd()  # returns to main canvas before defining pad2
    pad2 = TPad("pad2", "pad2", 0, 0.4, 1, 0.68)
    pad2.SetTopMargin(0)  # joins upper and lower plot
    pad2.SetGridx()
    pad2.SetBottomMargin(0)
    pad2.Draw()
    pad3 = TPad("pad3", "pad3", 0, 0.01, 1, 0.38)
    pad3.SetTopMargin(0)  # joins upper and lower plot
    pad3.SetGridx()
    pad3.SetBottomMargin(0)
    pad3.Draw()
    pad3.SetBottomMargin(0.4)
    pad1.cd()
    return c1, pad1, pad2, pad3


def get_histo_params(myvar):
    if 'base10logpt' in myvar:
        mybins = 10
        mymin = 1
        mymax = 4
    if 'bTagSF' in myvar:
        mybins = 20
        mymin = 0.8
        mymax = 1.5
    if 'X_resolved_WW_m' in myvar:
        mybins = 10
        mymin = 0
        mymax = 2000
    if 'nTrk' in myvar:
        mybins = 30
        mymin = 0
        mymax = 30
    if 'pdgid' in myvar:
        mybins = 25
        mymin = -2
        mymax = 23
    if 'truthFlavor' in myvar:
        mybins = 25
        mymin = -2
        mymax = 23
    if 'eta' in myvar:
        mybins = 20
        mymin = -3
        mymax = 3
    if 'qgtaggable' in myvar:
        mybins = 2
        mymin = 0
        mymax = 2
    if 'logpt' in myvar:
        mybins = 30
        mymin = 3
        mymax = 6.5
    if 'X_resolved_WW_m' in myvar:
        mybins = 10
        mymin = 0
        mymax = 2000
    if 'X_resolved_ZZ_m' in myvar:
        mybins = 10
        mymin = 0
        mymax = 2000
    if 'X_resolved_WZ_m' in myvar:
        mybins = 10
        mymin = 0
        mymax = 2000
    if 'pt' in myvar and 'log' not in myvar:
        mymin = 0
        if 'J1' in myvar:
            mybins = 30
            mymax = 600
        else:
            mymax = 300
            mybins = 15
    if 'tagJJ_m' in myvar:
        mybins = 20
        mymin = 0
        mymax = 5000
    if 'tagJJ_deta' in myvar:
        mybins = 10
        mymin = 0
        mymax = 10
    if 'pdgid_X_resolved_WW_m' == myvar:
        mybins = 10
        mymin = 0
        mymax = 2000

    return mybins, mymin, mymax

def get_max(histo_array):
    max_array = []
    for i in range(len(histo_array)):
        max_array.append(histo_array[i].GetMaximum())
    mymax = max(max_array)
    mymax = 1.2*mymax
    return mymax

def np_get_max(mymatrix):
    max_array=[]
    for i in range(mymatrix.shape[0]):
        for j in range(mymatrix.shape[1]):
            thismax = mymatrix[i][j].GetMaximum()
            max_array.append(thismax)
    mymax = max(max_array)
    return mymax

def get_var_title(var):
    if 'J1' in var:
        name = 'Leading Jet '
    if 'J2' in var:
        name = 'Sub-leading Jet '
    if 'pt' in var:
        return name + 'p_{T} [GeV]'
    if 'eta' in var:
        return name +'\eta'
    if 'nTrk' in var:
        return name + 'Number of Reconstructed Tracks'
    if 'pdgid' in var:
      return name + 'Truth Parton PDGID'
    if 'X_resolved' in var:
      return 'm_{lvqq} [GeV]'

def get_legend_hist_name(histname):
    if 'HVTWW_300' in histname:
        return 'Z\' 300GeV'
    if 'HVTWW_500' in histname:
        return 'Z\' 500GeV'
    if 'HVTWW_700' in histname:
        return 'Z\' 700GeV'
    if 'bkg' in histname:
        return 'Background'

def get_legend_values(myvar):
    if 'pt' in myvar or 'nTrk' in myvar:
        return 0.65,0.7,0.95,0.9
    if 'pdgid' in myvar:
        return 0.4, 0.7, 0.75, 0.9
    if 'eta' in myvar:
        return 0.78, 0.7, 0.98, 0.95

def make_2d_plots(myfile, sig, i, hfile, outputdir, j):
    #define selections
    if sig.eventqgselection[j] != "":
        selection = sig.weight + "*(" + sig.selection + ")*(" + sig.eventqgselection[j] + ")"
    else:
        selection = sig.weight + "*(" + sig.selection + ")"

    #obtain trees from file
    index1 = myfile.rfind("/") + 1
    index2 = myfile.rfind(".")
    filename = myfile[index1:index2]
    f = ROOT.TFile.Open(myfile, 'r')#add to only read file
    t = f.Get(sig.treename)

    tempf = ROOT.TFile.Open("temp.root", 'RECREATE')

    #create nTrk pt TH2F
    c3 = ROOT.TCanvas("c3", "c3", 0,  0, 1500, 1000)
    bins1, min1, max1 = get_histo_params(sig.twodvar1[i])
    bins2, min2, max2 = get_histo_params(sig.twodvar2[i])


    hist1 = ROOT.TH2F("hist1","", bins1, min1, max1, bins2, min2, max2)
    makehisto = sig.twodvar2[i] + ":" + sig.twodvar1[i] + ">>hist1"
    t.Draw(makehisto, selection, "goff")
    hist1.SetDirectory(tempf)

    #plotting aesthetics!
    canvasname = filename + "_" + sig.selection + "_" + sig.twodvar1[i] + "_" + sig.twodvar2[i] + "_" + sig.namesqgselection[j]
    hist1.SetTitle(canvasname)
    hist1.Draw("COLZ")

    axisname = sig.twodvar1[i]
    hist1.GetXaxis().SetTitle(sig.twodvar1[i])
    hist1.GetYaxis().SetTitle(sig.twodvar2[i])
    hist1.SetStats(0)
    c3.SetRightMargin(0.2)
    c3.Print(outputdir + "/" + canvasname + ".pdf")
    hfile.cd()
    c3.Write(canvasname)
    tempf.Write()
    tempf.Close()
    del hist1, c3, tempf
    return

def normalized_heat_maps_wrapper(sig, hfile, outputdir):
    cm = ROOT.TCanvas("cm", "cm", 0,  0, 1500, 1000)
    cm.SetRightMargin(0.2)
    for f in sig.files:
        for j in range(len(sig.ntrkvar)):
            if '300' in f:
                if j ==1:
                    hist300j1 = normalized_heat_maps_plots(f, j, sig, sig.heatmapvar, outputdir) 
                    hist300j1.Draw("COLZ")
                    cm.Print(outputdir + '/hvt300_sigj1.pdf')
            if '500' in f:
                if j ==1:
                    hist500j1 = normalized_heat_maps_plots(f, j, sig, sig.heatmapvar, outputdir) 
                    hist500j1.Draw("COLZ")
                    cm.Print(outputdir + '/hvt500_sigj1.pdf')
            if '700' in f:
                if j ==1:
                    hist700j1 = normalized_heat_maps_plots(f, j, sig, sig.heatmapvar, outputdir) 
                    hist700j1.Draw("COLZ")
                    cm.Print(outputdir + '/hvt700_sigj1.pdf')
            if 'bkg' in f:
                if j ==1:
                    bkgj1 = normalized_heat_maps_plots(f, j, sig, sig.heatmapvar, outputdir) 
                    bkgj1.Draw("COLZ")
                    cm.Print(outputdir + '/bkg_sigj1.pdf')

                    bkgj1.Add(hist300j1, -1)
                    bkgj1.Draw("COLZ")
                    f2 = ROOT.TF1("f2", "[0] + [1]*x + [2]*x + [3]*x", 0, 10)
                    f2.SetParameter(0, -16)
                    f2.SetParameter(1, 6)
                    f2.Draw("SAME")
                    cm.Print(outputdir + '/bkg_minus_signal_sigj1.pdf')


def normalized_heat_maps_plots(myfile, j, sig, var, outputdir):
    selection = sig.weight + "*(" + sig.selection + ")"#*(" + sig.taggablevar[j] + "==1)"
    #obtain trees from file
    index1 = myfile.rfind("/") + 1
    index2 = myfile.rfind(".")
    filename = myfile[index1:index2]
    f = ROOT.TFile.Open(myfile, 'r')#add to only read file
    t = f.Get(sig.treename)

    #create nTrk pt TH2F
    c3 = ROOT.TCanvas("c3", "c3", 0,  0, 1500, 1000)
    logptbins, logptmin, logptmax = get_histo_params('logpt')
    #logptbins, logptmin, logptmax = get_histo_params('base10logpt')
    ntrkbins, ntrkmin, ntrkmax = get_histo_params('nTrk')
    if var=="pt" or var=="e":
        variables = sig.ntrkvar[j] + ":TMath::Log(" + sig.ptvar[j] +  ")"

    hist1 = ROOT.TH2F("hist1","", logptbins, logptmin, logptmax, ntrkbins, ntrkmin, ntrkmax)
    makehisto = variables + ">>hist1"
    t.Draw(makehisto, selection, "goff")
    hist1.SetDirectory(0)

    canvasname = filename + "_" + sig.selection + "_" + sig.ntrkvar[j]
    hist1.SetTitle(canvasname)

    if var =="e" or var == "pt":
        axisname = "ln(" + var + ")"
    if var == "eta":
        axisname ="eta"
    hist1.GetXaxis().SetTitle(axisname)
    hist1.GetYaxis().SetTitle("nTrk")
    hist1.SetStats(0)
    hist1.Scale(1/hist1.GetSumOfWeights())

    return hist1

def make_quark_ratio(myfile, sig, hfile, outputdir, j):
    if sig.eventqgselection[j] != "":
        selection = sig.weight + "*(" + sig.selection + ")*(" + sig.eventqgselection[j] + ")"
    else:
        selection = sig.weight + "*(" + sig.selection + ")"

    #obtain trees from file
    index1 = myfile.rfind("/") + 1
    index2 = myfile.rfind(".")
    filename = myfile[index1:index2]
    f = ROOT.TFile.Open(myfile, 'r')#add to only read file
    t = f.Get(sig.treename)

    tempf = ROOT.TFile.Open("temp.root", 'RECREATE')

    #create nTrk pt TH2F
    c3 = ROOT.TCanvas("c3", "c3", 0,  0, 1500, 1000)
    bins1, min1, max1 = get_histo_params(sig.eventvar)

if __name__ == '__main__': 
    yieldfunctions()

