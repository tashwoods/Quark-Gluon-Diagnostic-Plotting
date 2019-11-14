class signalObject:
    def __init__(self, signal, slope, inter):
        self.signal = signal
        self.weight = 'weight'
        self.treename = 'nominal'
        self.qgslope = slope
        self.qginter = inter
        self.ntrkmin = '0'

        #Set directories that store input files
        if '1lep' in signal:
            self.basedir = '/gpfs/slac/atlas/fs1/d/woodsn/data/qgthesis_files/files/'
        else:
            self.basedir = '/gpfs/slac/atlas/fs1/d/woodsn/StatFrameWork/VVSemileptonicStats/SlimmedNtuples/2lep_jan31_all/Nominal/data/itest01/'


        #Setup plotting variables and their constraints
        if signal == '1lep_HVTWW':
            self.selection = 'Pass_Res_GGF_WW_SR'
            files = ['HVTWW_300_1lep.root', 'HVTWW_500_1lep.root', 'allbkg-ade.root'] #'HVTWW_700_1lep.root'] #'allbkg-ade.root']
            self.shortfilenames=['Z\' 300GeV', 'Z\' 500GeV', 'Backgrounds']
            self.indices = [0,1,2]
            self.masses = ['300', '500', 'na']
            self.color = [8,63,53]
            self.filesxs = [3.172,4.6140E-01,1.0]
            self.eventvar = 'X_resolved_WW_m'
            self.binning = [0, 260, 300, 330, 370, 410, 450, 490, 540, 590, 640, 700, 760, 820, 890, 960, 1040]
            self.variables = ['sigWJ1_nTrk', 'sigWJ1_pdgid', 'sigWJ1_pt', 'sigWJ1_eta']
            self.jets = ['sigWJ1']
            self.ntrkvar = ['sigWJ1_nTrk']
            self.ptvar = ['sigWJ1_pt']
            self.truthvar = ['sigWJ1_pdgid']
            self.taggablevar = ['1']#'sigWJ1_qgtaggable', 'sigWJ2_qgtaggable']
            self.heatmapvar = 'pt'
            self.twodvar1 = ['X_resolved_WW_m']
            self.twodvar2 = ['sigWJ1_pdgid']
        elif signal == '1lep_HVTWZ':
            self.selection = 'Pass_Res_GGF_WZ_Tag_SR'
            files = ['HVTWZ_300_SR.root', 'HVTWZ_700_SR.root', 'bkg_SR.root']
            self.filesxs = [13.8, 10, 1]
            self.eventvar = 'X_resolved_WW_m'
            self.binning = [0, 300, 360, 420, 500, 575, 660, 755, 860, 975, 1100, 1500, 2000]
            self.variables = ['sigZJ1_nTrk', 'sigZJ2_nTrk', 'sigZJ1_pdgid', 'sigZJ2_pdgid', 'sigZJ1_eta', 'sigZJ2_eta', 'sigZJ1_qgtaggable', 'sigZJ2_qgtaggable', 'X_resolved_WW_m', 'sigZJ1_pt', 'sigZJ2_pt', 'X_resolved_WW_m', 'weightBtagSF']
            self.jets = ['sigZJ1', 'sigZJ2']
            self.ntrkvar = ['sigZJ1_nTrk', 'sigZJ2_nTrk']
            self.ptvar = ['sigZJ1_pt', 'sigZJ2_pt']
            self.truthvar = ['sigZJ1_pdgid', 'sigZJ2_pdgid']
            self.taggablevar = ['sigZJ1_qgtaggable', 'sigZJ2_qgtaggable']
            self.heatmapvar = 'pt'
            self.twodvar1 = ['tagJJ_m']
            self.twodvar2 = ['tagJJ_deta']
        else:
            self.twodvar1 = [ 'tagJJ_m']
            self.twodvar2 = ['tagJJ_deta']
            self.selection = 'Pass_Res_VBF_WW_SR'
            self.eventvar = 'X_resolved_WW_m'
            files = ['HVTWWVBF_300_SR.root', 'HVTWWVBF_700_SR.root', 'bkg_SR.root']
            self.filesxs = [20, 280, 1]
            self.jets = ['sigWJ1', 'sigWJ2', 'tagJ1', 'tagJ2']
            self.ntrkvar = ['sigWJ1_nTrk', 'sigWJ2_nTrk', 'tagJ1_nTrk', 'tagJ2_nTrk']
            self.ptvar = ['sigWJ1_pt', 'sigWJ2_pt', 'tagJ1_pt', 'tagJ2_pt']
            self.truthvar = ['sigWJ1_pdgid', 'sigWJ2_pdgid', 'tagJ1_pdgid', 'tagJ2_pdgid']
            self.taggablevar = ['sigWJ1_qgtaggable', 'sigWJ2_qgtaggable', 'tagJ1_qgtaggable', 'tagJ2_qgtaggable']
            self.eventvar = 'X_resolved_WW_m'
            self.binning = [0, 300, 360, 420, 500, 575, 660, 755, 860, 975, 1100, 1500, 2000]
            self.variables = ['sigWJ1_nTrk', 'sigWJ2_nTrk', 'sigWJ1_pdgid', 'sigWJ2_pdgid', 'tagJ1_nTrk', 'tagJ2_nTrk', 'tagJ1_pdgid', 'tagJ2_pdgid', 'sigWJ1_qgtaggable', 'sigWJ2_qgtaggable', 'tagJ1_qgtaggable', 'tagJ2_qgtaggable', 'X_resolved_WW_m', 'sigWJ1_eta', 'sigWJ2_eta', 'tagJ1_eta', 'tagJ2_eta', 'sigWJ1_pt', 'sigWJ2_pt', 'tagJ1_pt', 'tagJ2_pt', 'tagJJ_m', 'tagJJ_deta']
            self.heatmapvar = 'pt'

        #Append basedir to filesnames
        self.files = [ self.basedir + s for s in files ]

        if 'VBF' in signal:
            self.eventqgselection = ['', 'Pass_Res_qg_tag_sigjets', '(Pass_Res_qg_tag_sigjets && Pass_Res_qg_tag_tagjets)']
            self.namesqgselection = ['Nominal', 'tag_sigjets', 'tag_sigjets_tagjets']
        else:
            self.qgptmin = '50'
            sigWJ1_taggable = "((sigWJ1_pt > " + self.qgptmin + ") && (fabs(sigWJ1_eta)<2.1))"
            sigWJ2_taggable = "((sigWJ2_pt > " + self.qgptmin + ") && (fabs(sigWJ2_eta)<2.1))"

            self.nottaggable1 = '(!' + sigWJ1_taggable + ')'
            self.taggable_tagged1 = '(' + sigWJ1_taggable + ' && sigWJ1_nTrk > ' + str(self.ntrkmin) + ' && sigWJ1_nTrk < (' + str(self.qgslope) + '*TMath::Log(sigWJ1_pt)+' + str(self.qginter) + '))'
            self.sigJ1_rqmt = '(' + self.nottaggable1 + ' || ' + self.taggable_tagged1 + ')'

            self.nottaggable2 = '((sigWJ2_qgtaggable == 0) || (sigWJ2_qgtaggable==1 && sigWJ2_pt < ' + str(self.qgptmin) + '))'
            self.taggable_tagged2 = '(sigWJ2_qgtaggable ==1 && sigWJ2_pt > ' + str(self.qgptmin) + ' && sigWJ2_nTrk > ' + str(self.ntrkmin) + ' && sigWJ2_nTrk < (' + str(self.qgslope) + '*TMath::Log(sigWJ2_pt)+' + str(self.qginter) + '))'
            self.sigJ2_rqmt = '(' + self.nottaggable2 + ' || ' + self.taggable_tagged2 + ')'

            thisselection = '( ' + self.sigJ1_rqmt + ')' #+ ' && ' + self.sigJ2_rqmt + ')'
            #thisselection = '((sigWJ1_qgtaggable==1 && sigWJ1_nTrk > 0 && (sigWJ1_nTrk < (' + str(self.qgslope) + '*TMath::Log(sigWJ1_pt)+' + str(self.qginter) + ')) || sigWJ1_qgtaggable == 0 || (sigWJ1_qgtaggable==1 && sigWJ1_pt < ' + str(self.qgptmin) + ')) && (sigWJ2_qgtaggable==1 && sigWJ2_nTrk > 0 && (sigWJ2_nTrk < (' + str(self.qgslope) + '*TMath::Log(sigWJ2_pt)+ ' + str(self.qginter) + ')) || sigWJ2_qgtaggable==0 || (sigWJ2_qgtaggable==1 && sigWJ2_pt < ' + str(self.qgptmin) + ')) '
            self.eventqgselection = ['', thisselection]
            self.namesqgselection = ['Nominal', 'tag_sigjets']
