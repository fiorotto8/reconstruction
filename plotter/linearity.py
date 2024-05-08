import os, math, optparse, ROOT
ROOT.gROOT.SetBatch(1)
from simple_plot import getCanvas, doLegend
from mcPlots import doTinyCmsPrelim

ROOT.gROOT.LoadMacro('./fitter/libCpp/RooDoubleCBFast.cc+')
ROOT.gROOT.LoadMacro('./fitter/libCpp/RooCruijff.cc+')

def full_error_e(energy,staterr,deltax=2):
    # from the z-scan, we observe ~10% response change / 5cm for z ~ 21 cm
    syst_pos_rel = 0.02 * deltax
    syst_pos = syst_pos_rel * energy
    return math.hypot(syst_pos,staterr)

def ic_delta_multisource(mat,k=1.2):
    return k if mat!='Fe' else 1

def fitOne(filein,selection,erange,mrange,suffix,energies,nbins,pdir):

    var = 'sc_integral*6/8000'
    # fitting Fe peak in Ti transparency mode
    if suffix=='Ti':
        var="{base} * {calib}".format(base=var,calib=5.9/6.15)
    
    xmin,xmax=erange
    mark = ROOT.kFullSquare
    line = ROOT.kAzure-6

    rf = ROOT.TFile.Open(filein)
    t = rf.Get('Events')
    spectrum = ROOT.TH1F('spectrum','',nbins,xmin,xmax)
    spectrum.SetMarkerStyle(mark)
    spectrum.SetLineColor(ROOT.kBlack)
    t.Draw('{var}>>spectrum'.format(var=var),'({var}>{xmin}) * ({var}<{xmax}) * ({sel})'.format(var=var,xmin=xmin,xmax=xmax,sel=selection))
    deltae = energies[1]-energies[0]
    av=energies[0]

    work = ROOT.RooWorkspace()
    work.factory("mean1[{av},{mmin},{mmax}]".format(av=av,mmin=mrange[0],mmax=mrange[1]))
    if suffix in ['XXX']:
        work.factory('CBShape::cb1(x[{xmin},{xmax}],mean1,sigma[1,0.5,3],alpha[1,0.1,10],n[9,1,20])'.format(xmin=xmin,xmax=xmax))
    elif suffix in ['Cu','Fe','Rb','Mo','Ag','Ba','Ti','Ca']:
        work.factory('Cruijff::cb1(x[{xmin},{xmax}],mean1,sigmaL[1,0.5,4],sigma[1,0.2,4],alphaL[1,0.1,10],alphaR[0,0.1,10])'.format(xmin=xmin,xmax=xmax))
    else:
        work.factory('DoubleCBFast::cb1(x[{xmin},{xmax}],mean1,sigma[1,0.8,10],alpha1[1,0.01,10],n1[5,1,20],alpha2[1,0.01,10],n2[5,1,20])'.format(xmin=xmin,xmax=xmax))
    work.factory("expr::mean2('mean1+delta',mean1,delta[{delta}])".format(delta=deltae))
    if suffix in ['XXX']:
        work.factory('CBShape::cb2(x[{xmin},{xmax}],mean2,sigma,alpha,n)'.format(xmin=xmin,xmax=xmax))
    elif suffix in ['Cu','Fe','Rb','Mo','Ag','Ba','Ti','Ca']:
        work.factory('Cruijff::cb2(x[{xmin},{xmax}],mean2,sigmaL,sigma,alphaL,alphaR)'.format(xmin=xmin,xmax=xmax))
    else:
        work.factory('DoubleCBFast::cb2(x[{xmin},{xmax}],mean2,sigma,alpha1,n1,alpha2,n2)'.format(xmin=xmin,xmax=xmax))
    if suffix in ['Ti','Ca']:
        # this is to model the Fe55 peak
        work.factory('Cruijff::cb3(x[{xmin},{xmax}],mean3[5.9,5,8],sigmaL,sigma,alphaL,alphaR)'.format(xmin=xmin,xmax=xmax))
        work.factory("expr::mean4('mean3+deltaFe',mean3,deltaFe[0.6])")
        work.factory('Cruijff::cb4(x[{xmin},{xmax}],mean4,sigmaL,sigma,alphaL,alphaR)'.format(xmin=xmin,xmax=xmax))


    work.factory("a[0,-100,1000]")
    work.factory("b[0,-100,100]")
    work.factory("c[0,-100,100]")
    work.factory("d[0,-100,100]")
    work.factory("e[0,-100,100]")
    if suffix in ['Fe']:
        work.factory('Bernstein::bkg(x,{a})')
    else:
        work.factory('Bernstein::bkg(x,{a,b,c,d,e})')
    
    x = work.var('x')

    work.factory('nSig1[{ns},0,1e9]'.format(ns=0.2*spectrum.Integral()))
    work.factory("expr::nSig2('nSig1*frac',nSig1,frac[0.1,0,0.3])")
    work.factory('nBkg[{ns}]'.format(ns=0.5*spectrum.Integral()))
    if suffix in ['Ti','Ca']:
        work.factory('nSig3[{ns},0,1e9]'.format(ns=0.7*spectrum.Integral()))
        work.factory("expr::nSig4('nSig3*frac2',nSig3,frac2[0.05,0.0,0.4])")
        work.factory("SUM::pdfTot(nSig1*cb1,nSig2*cb2,nSig3*cb3,nSig4*cb4,nBkg*bkg)");
    else:
        work.factory("SUM::pdfTot(nSig1*cb1,nSig2*cb2,nBkg*bkg)");

    if suffix=='Fe':
        work.var('nBkg').setVal(0)
        work.var('nBkg').setConstant(ROOT.kTRUE)
    if suffix=='Ba':
        work.var('sigma').setRange(1,6)
    elif suffix=='Cu':
        work.var('sigma').setRange(0.5,1)
    if suffix in ['Ba']:
        work.var('alphaR').setRange(0.,0.3)
    if suffix in ['Rb']:
         work.var('alphaR').setVal(1)
    elif suffix in ['Cu']:
        work.var('alphaR').setRange(0.1,10)
        work.var('alphaR').setConstant(ROOT.kFALSE)
    elif suffix in ['Ag']:
        work.var('alphaR').setVal(0.1)
        work.var('alphaR').setRange(0.1,0.4)
        work.var('sigma').setRange(1.8,3.1)
    if suffix=='Ti':
        work.var('frac2').setVal(0.1);       work.var('frac2').setConstant(ROOT.kTRUE)
        work.var('frac').setRange(0.3,0.5);  work.var('frac').setConstant(ROOT.kTRUE)
        work.var('alphaL').setVal(0.251);    #work.var('alphaL').setConstant(ROOT.kTRUE)
        work.var('alphaR').setVal(0.252);    #work.var('alphaR').setConstant(ROOT.kTRUE)
        work.var('sigma').setVal(0.616);     #work.var('sigma').setConstant(ROOT.kTRUE)
        work.var('sigmaL').setRange(0.2,1.0);    #work.var('sigmaL').setConstant(ROOT.kTRUE)
        work.var('mean3').setVal(5.33);      #work.var('mean3').setConstant(ROOT.kTRUE)
        work.var('mean1').setVal(4.5)
        
        work.var('a').setVal(2.62289e-02);  work.var('a').setConstant(ROOT.kTRUE)
        work.var('b').setVal(8.81561e-02);  work.var('b').setConstant(ROOT.kTRUE)
        work.var('c').setVal(6.68959e-04);  work.var('c').setConstant(ROOT.kTRUE)
        work.var('d').setVal(2.18863e-02);  work.var('d').setConstant(ROOT.kTRUE)
        work.var('e').setVal(2.00491e-02);  work.var('e').setConstant(ROOT.kTRUE)
    if suffix=='Ca':
        work.var('frac2').setRange(0.2,0.3); work.var('frac').setRange(0.1,0.3); 
        work.var('alphaL').setVal(0.251);    work.var('alphaL').setConstant(ROOT.kTRUE)
        work.var('alphaR').setVal(0.);       #work.var('alphaR').setConstant(ROOT.kTRUE)
        work.var('sigma').setVal(0.616);     #work.var('sigma').setRange(0.4,1.3)
        work.var('sigmaL').setVal(0.616);    #work.var('sigmaL').setRange(0.4,1.3)
        work.var('mean3').setVal(5.33);      work.var('mean3').setConstant(ROOT.kTRUE)
        work.var('mean1').setVal(3.5)
        
        work.var('a').setVal(2.62289e-02);  work.var('a').setConstant(ROOT.kTRUE)
        work.var('b').setVal(8.81561e-02);  work.var('b').setConstant(ROOT.kTRUE)
        work.var('c').setVal(6.68959e-04);  work.var('c').setConstant(ROOT.kTRUE)
        work.var('d').setVal(2.18863e-02);  work.var('d').setConstant(ROOT.kTRUE)
        work.var('e').setVal(2.00491e-02);  work.var('e').setConstant(ROOT.kTRUE)


        
    work.Print()
    
    rooData = ROOT.RooDataHist("spectrum","spectrum",ROOT.RooArgList(work.var("x")),spectrum)
    getattr(work,'import')(rooData)

    frame = x.frame()
    frame.SetTitle('')
    
    # fit spectrum
    rooData.plotOn(frame,ROOT.RooFit.MarkerStyle(mark),ROOT.RooFit.MarkerSize(2),ROOT.RooFit.Name('data'))
    pdf = work.pdf('pdfTot')
    pdf.fitTo(rooData,ROOT.RooFit.Save(),ROOT.RooFit.PrintLevel(10))
    pdf.plotOn(frame,ROOT.RooFit.LineColor(line), ROOT.RooFit.Name('ptot'))
    pdf.plotOn(frame,ROOT.RooFit.Components("bkg"),ROOT.RooFit.LineStyle(ROOT.kDotted), ROOT.RooFit.LineColor(ROOT.kCyan), ROOT.RooFit.Name('pbkg'))
    pdf.plotOn(frame,ROOT.RooFit.Components("cb1"), ROOT.RooFit.LineStyle(ROOT.kDotted), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name('pcb1'))
    pdf.plotOn(frame,ROOT.RooFit.Components("cb2"), ROOT.RooFit.LineStyle(ROOT.kDotted), ROOT.RooFit.LineColor(ROOT.kOrange), ROOT.RooFit.Name('pcb2'))
    if suffix in ['Ti','Ca']:
        pdf.plotOn(frame,ROOT.RooFit.Components("cb3"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kSpring+4), ROOT.RooFit.Name('pcb3'))
        pdf.plotOn(frame,ROOT.RooFit.Components("cb4"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kSpring+1), ROOT.RooFit.Name('pcb4'))        
    rooData.plotOn(frame,ROOT.RooFit.MarkerStyle(mark))

    frame.GetYaxis().SetTitle("superclusters")
    frame.GetXaxis().SetTitle("Energy (keV)")
    frame.GetXaxis().SetTitleOffset(1.2)
    frame.GetXaxis().SetTitleFont(42)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetTitleOffset(1.1)
    frame.GetXaxis().SetLabelFont(42)
    frame.GetXaxis().SetLabelSize(0.05)
    frame.GetXaxis().SetLabelOffset(0.007)
    frame.GetYaxis().SetTitleFont(42)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetTitleOffset(1.3)
    frame.GetYaxis().SetLabelFont(42)
    frame.GetYaxis().SetLabelSize(0.05)
    frame.GetYaxis().SetLabelOffset(0.007)
    frame.GetXaxis().SetNdivisions(510)

    c = getCanvas()
    frame.Draw()

    # Add legend
    leg = ROOT.TLegend(0.6,0.6,0.85,0.88)
    leg.SetFillStyle(0)
    leg.SetLineColor(0)
    leg.SetTextSize(0.02)
    leg.SetTextFont(42)
    leg.AddEntry("pdata","data with {source}".format(source=suffix),"ep")
    leg.AddEntry("tot","Parametric Total Model","L")
    leg.AddEntry("pcb1","Signal Model 1","L")
    leg.AddEntry("pcb2","Signal Model 2","L")
    if suffix in ['Ti','Ca']:
        leg.AddEntry("pcb3","^{55}Fe Signal Model 1","L")
        leg.AddEntry("pcb4","^{55}Fe Signal Model 2","L")        
    leg.AddEntry("pbkg","Background Model","L")
    leg.Draw("Same")

    doTinyCmsPrelim("#bf{CYGNO}","(LIME)",lumi=0,textSize=0.05,xoffs=-0.06)
    
    for ext in ['pdf','png','root']:
        c.SaveAs('{pdir}/energy_{suff}.{ext}'.format(pdir=pdir,ext=ext,suff=suffix))

    m    = work.var('mean1').getVal()
    merr = work.var('mean1').getError()
    s    = work.var('sigma').getVal()
    serr = work.var('sigma').getError()
    if suffix in ['Ti','Ca']:
        m2 = work.var('mean3').getVal()
        m2err = work.var('mean3').getError()
    else:
        m2 = None
        m2err = None

    return [m,merr,s,serr,m2,m2err]

if __name__ == "__main__":

    parser = optparse.OptionParser(usage='usage: %prog [opts] ', version='%prog 1.0')
    parser.add_option("--pdir", "--print-dir", dest="printDir", type="string", default="plots", help="print out plots in this directory");
    (options, args) = parser.parse_args()

    outname = options.printDir
    if not os.path.exists(outname):
        os.system("mkdir -p "+outname)
        if os.path.exists("../index.php"): os.system("cp ../index.php "+outname)
    print("Will save plots to ",outname)
    
    base_sel = "(TMath::Hypot(sc_xmean-2304/2,sc_ymean-2304/2)<800)"

    materials = ["Ag","Ba","Cu","Rb","Mo","Fe","Ti","Ca"]

    energies = {"Ca" : [3.69,  4.01],
                "Ti" : [4.51,  4.93],
                "Fe" : [5.9,   6.5],
                "Cu" : [8.04,  8.91],
                "Rb" : [13.37, 14.97],
                "Mo" : [17.44, 19.63],
                "Ag" : [22.10, 24.99],
                "Ba" : [32.06, 36.55],
                }
        
    
    extras = {"Ca" : "(0.152*sc_length < 80)  * (sc_rms>7.5) * (run!=6153 && run!=6164 && run!=6175 && run!=6186 && run!=6197 && run!=6208 && run!=6221 && run!=6232 && run!=6243 && run!=6254 && run!=6265 && run!=6276 && run!=6287 && run!=6290)",
              "Fe" : "(0.152*sc_length < 80)  * (sc_rms>7.5)",
              "Cu" : "(0.152*sc_length < 80)  * (sc_rms>7.5) * (run>=5801 && run<=5805)",
              "Rb" : "(0.152*sc_length < 80)  * (sc_rms>6.5) * (run>=5806 && run<=5810)",
              "Mo" : "(0.152*sc_length < 80)  * (sc_rms>6.5) * (run>=5811 && run<=5820)",
              "Ag" : "(0.152*sc_length < 80)  * (sc_rms>6.5) * (run>=5821 && run<=5830)",
              "Ba" : "(0.152*sc_length < 120) * (sc_rms>6.5) * (run>=5831 && run<=5845)",
              "Ti" : "(0.152*sc_length < 80)  * (sc_rms>7.5) * (run>=6122 && run<=6131)",
              }


    eranges = {"Fe" : (2,12),
               "Cu" : (2,16),
               "Rb" : (1,25),
               "Mo" : (8,30),
               "Ag" : (9,40),
               "Ba" : (11,60),
               "Ti" : (1,11),
               "Ca" : (0.5,16),}

    nbins = {"Fe" : 120,
             "Cu" : 45,
             "Rb" : 45,
             "Mo" : 55,
             "Ag" : 50,
             "Ba" : 50,
             "Ti" : 60,
             "Ca" : 70}

    mrange = {"Fe" : (5,8),
              "Cu" : (6,9),
              "Rb" : (9,11.5),
              "Mo" : (12,18),
              "Ag" : (17,22),
              "Ba" : (25,36),
              "Ti" : (3.8,5.0),
              "Ca" : (3,4.5),}

    rfiles = {"Fe" : "reco_fe_26cm.root",
              "Cu" : "reco_multisource.root",
              "Rb" : "reco_multisource.root",
              "Mo" : "reco_multisource.root",
              "Ag" : "reco_multisource.root",
              "Ba" : "reco_multisource.root",
              "Ti" : "reco_titanium.root",
              "Ca" : "reco_calcium.root",}
              
    respdic = {}
    respdic["Fe"] = [energies['Fe'][0],0,0.9,0]

    for mat in materials:
        respdic[mat] = fitOne("trees-lnf/%s" % rfiles[mat],'({base})*({extra})'.format(base=base_sel,extra=extras[mat]),eranges[mat],mrange[mat],mat,energies[mat],nbins[mat],outname)

    calibFe = energies['Fe'][0]/respdic['Fe'][0]

    
    fout = ROOT.TFile.Open("linearity.root","recreate")
    resp  = ROOT.TGraphErrors(len(energies))
    resp2 = ROOT.TGraphErrors(len(energies))
    reso = ROOT.TGraphErrors(len(energies))
    i=0
    for mat,etrue in energies.items():
        response1 =  respdic[mat][0]*calibFe*ic_delta_multisource(mat)
        response2 = (respdic[mat][0]+energies[mat][1]-energies[mat][0])*calibFe*ic_delta_multisource(mat)
        unc1      = full_error_e(response1,respdic[mat][1]*calibFe*ic_delta_multisource(mat))
        unc2      = full_error_e(response2,respdic[mat][1]*calibFe*ic_delta_multisource(mat))
        if mat in ['Ti','Ca']: # in these cases the Fe peak dominates the fit for the sigmaR
            scaleForResolNorm = respdic[mat][4]*calibFe*ic_delta_multisource(mat)
        else:
            scaleForResolNorm = respdic[mat][0]*calibFe*ic_delta_multisource(mat)
            
        resp.SetPoint(i,energies[mat][0],response1)
        resp.SetPointError(i,0,unc1)
        resp2.SetPoint(i,energies[mat][1],response2)
        resp2.SetPointError(i,0,unc2)
        reso.SetPoint(i,energies[mat][0],100*respdic[mat][2]/scaleForResolNorm/ic_delta_multisource(mat))
        reso.SetPointError(i,0,100*respdic[mat][3]/scaleForResolNorm/ic_delta_multisource(mat))
        i+=1

    
    ## response linearity
    c = getCanvas()
    resp.SetTitle('')
    resp.Draw("AP")
    resp.SetMarkerStyle(ROOT.kFullCircle)
    resp2.SetMarkerColor(ROOT.kAzure-3)
    resp2.SetLineColor(ROOT.kAzure-3)
    resp.SetMarkerSize(2)
    resp.GetXaxis().SetLimits(0,50)
    resp.GetXaxis().SetTitle("Expected energy (keV)")
    resp.GetYaxis().SetLimits(0,8e4)
    resp.GetYaxis().SetRangeUser(0,50)
    resp.GetYaxis().SetTitle("Measured energy (keV)")
    resp.GetXaxis().SetTitleOffset(1.2)
    resp.GetXaxis().SetTitleFont(42)
    resp.GetXaxis().SetTitleSize(0.05)
    resp.GetXaxis().SetTitleOffset(1.1)
    resp.GetXaxis().SetLabelFont(42)
    resp.GetXaxis().SetLabelSize(0.05)
    resp.GetXaxis().SetLabelOffset(0.007)
    resp.GetYaxis().SetTitleFont(42)
    resp.GetYaxis().SetTitleSize(0.05)
    resp.GetYaxis().SetTitleOffset(1.2)
    resp.GetYaxis().SetLabelFont(42)
    resp.GetYaxis().SetLabelSize(0.05)
    resp.GetYaxis().SetLabelOffset(0.007)
    resp.GetXaxis().SetNdivisions(510)

    
    resp2.Draw("P")
    resp2.SetMarkerStyle(ROOT.kFullSquare)
    resp2.SetMarkerColor(ROOT.kOrange+7)
    resp2.SetLineColor(ROOT.kOrange+7)
    resp2.SetMarkerSize(2)

    # Add legend
    leg = ROOT.TLegend(0.2,0.75,0.4,0.85)
    leg.SetTextFont(42)
    leg.SetFillStyle(0)
    leg.SetLineColor(0)
    leg.SetTextSize(0.04)
    leg.AddEntry(resp,"1^{st} line","ep")
    leg.AddEntry(resp2,"2^{nd} line","ep")
    leg.Draw("Same")

    #resp.Fit("pol1")

    line = ROOT.TLine(0,0,50,50)
    line.SetLineColor(ROOT.kBlue)
    line.SetLineStyle(ROOT.kDashed)
    line.Draw()

    doTinyCmsPrelim("#bf{CYGNO}","(LIME)",lumi=0,textSize=0.05,xoffs=-0.06)

    fout.cd()
    resp.Write()
    c.Write()
    for ext in ['pdf','png','root']:
        c.SaveAs('{pdir}/linearity.{ext}'.format(pdir=outname,ext=ext))
    fout.Close()
    
    # energy resolution
    c = getCanvas()
    reso.SetTitle('')
    reso.Draw("AP")
    reso.SetMarkerStyle(ROOT.kFullCircle)
    reso.SetMarkerSize(2)
    reso.GetXaxis().SetLimits(0,50)
    reso.GetXaxis().SetTitle("Energy (keV)")
    reso.GetYaxis().SetLimits(0,30)
    reso.GetYaxis().SetRangeUser(0,30)
    reso.GetYaxis().SetTitle("#sigma_{E}/E (%)")        
    doTinyCmsPrelim("#bf{CYGNO}","(LIME)",lumi=0,textSize=0.04,xoffs=-0.05)
    for ext in ['pdf','png','root']:
        c.SaveAs('{pdir}/resolution.{ext}'.format(pdir=outname,ext=ext))




