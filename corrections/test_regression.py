# example USAGE: python test_regression.py ../prod-winter23-20230214/Merged/merged_feruns_8882_9857.root params_gbrtrain.txt --loadPanda regrdata_lngs_run2.pkl
import ROOT
ROOT.gROOT.SetBatch(True)

import numpy as np
import math, os, joblib

from gbr_trainer import GBRLikelihoodTrainer,fill_hist
ROOT.gROOT.LoadMacro('../plotter/fitter/libCpp/RooCruijff.cc+')

Zsteps = [5,15,25,36,48]

def doSpam(text,x1,y1,x2,y2,align=12,fill=False,textSize=0.033,_noDelete={}):
    cmsprel = ROOT.TPaveText(x1,y1,x2,y2,"NDC");
    cmsprel.SetTextSize(textSize);
    cmsprel.SetFillColor(0);
    cmsprel.SetFillStyle(1001 if fill else 0);
    cmsprel.SetLineStyle(2);
    cmsprel.SetLineColor(0);
    cmsprel.SetTextAlign(align);
    cmsprel.SetTextFont(42);
    cmsprel.AddText(text);
    cmsprel.Draw("same");
    _noDelete[text] = cmsprel; ## so it doesn't get deleted by PyROOT
    return cmsprel

def doTinyCmsPrelim(textLeft="#bf{CYGNO}",textRight="(LIME - Run2)",hasExpo=False,textSize=0.033,lumi=None, xoffs=0, options=None, doWide=False):
    if textLeft  == "_default_": textLeft  = "#bf{CYGNO}"
    if textRight == "_default_": textRight = "(LIME)"
    if lumi      == None       : lumi      = 1
    if   lumi > 3.54e+1: lumitext = "%.0f fb^{-1}" % lumi
    elif lumi > 3.54e+0: lumitext = "%.1f fb^{-1}" % lumi
    elif lumi > 3.54e-1: lumitext = "%.2f fb^{-1}" % lumi
    elif lumi > 3.54e-2: lumitext = "%.0f pb^{-1}" % (lumi*1000)
    elif lumi > 3.54e-3: lumitext = "%.1f pb^{-1}" % (lumi*1000)
    else               : lumitext = "%.2f pb^{-1}" % (lumi*1000)
    lumitext = "%.1f fb^{-1}" % lumi
    textLeft = textLeft.replace("%(lumi)",lumitext)
    textRight = textRight.replace("%(lumi)",lumitext)
    if textLeft not in ['', None]:
        doSpam(textLeft, (.28 if hasExpo else 0.2 if doWide else .12)+xoffs, .94, .60+xoffs, .94, align=12, textSize=textSize)
    if textRight not in ['', None]:
        doSpam(textRight,(0.5 if doWide else .55)+xoffs, .94, .82+xoffs if doWide else .91+xoffs, .94, align=32, textSize=textSize)

def getCanvas(name='c'):

    ROOT.gStyle.SetPalette(ROOT.kRainBow)
    ROOT.gStyle.SetNumberContours(51)
    ROOT.gErrorIgnoreLevel = 100
    ROOT.gStyle.SetOptStat(0)

    c = ROOT.TCanvas(name,'',1200,1200)
    lMargin = 0.14
    rMargin = 0.10
    bMargin = 0.15
    tMargin = 0.10
    c.SetLeftMargin(lMargin)
    c.SetRightMargin(rMargin)
    c.SetTopMargin(tMargin)
    c.SetBottomMargin(bMargin)
    c.SetFrameBorderMode(0);
    c.SetBorderMode(0);
    c.SetBorderSize(0);
    return c

def doLegend(histos,labels,styles,corner="TR",textSize=0.035,legWidth=0.18,legBorder=False,nColumns=1):
    nentries = len(histos)
    (x1,y1,x2,y2) = (.85-legWidth, .7 - textSize*max(nentries-3,0), .90, .89)
    if corner == "TR":
        (x1,y1,x2,y2) = (.85-legWidth, .7 - textSize*max(nentries-3,0), .90, .89)
    elif corner == "TC":
        (x1,y1,x2,y2) = (.5, .7 - textSize*max(nentries-3,0), .5+legWidth, .89)
    elif corner == "TL":
        (x1,y1,x2,y2) = (.2, .7 - textSize*max(nentries-3,0), .2+legWidth, .89)
    elif corner == "BR":
        (x1,y1,x2,y2) = (.85-legWidth, .15 + textSize*max(nentries-3,0), .90, .25)
    elif corner == "BC":
        (x1,y1,x2,y2) = (.5, .15 + textSize*max(nentries-3,0), .5+legWidth, .25)
    elif corner == "BL":
        (x1,y1,x2,y2) = (.2, .23 + textSize*max(nentries-3,0), .33+legWidth, .35)
    leg = ROOT.TLegend(x1,y1,x2,y2)
    leg.SetNColumns(nColumns)
    leg.SetFillColor(0)
    leg.SetFillColorAlpha(0,0.6)  # should make the legend semitransparent (second number is 0 for fully transparent, 1 for full opaque)
    #leg.SetFillStyle(0) # transparent legend, so it will not cover plots (markers of legend entries will cover it unless one changes the histogram FillStyle, but this has other effects on color, so better not touching the FillStyle)
    leg.SetShadowColor(0)
    if not legBorder:
        leg.SetLineColor(0)
        leg.SetBorderSize(0)  # remove border  (otherwise it is drawn with a white line, visible if it overlaps with plots
    leg.SetTextFont(42)
    leg.SetTextSize(textSize)
    for (plot,label,style) in zip(histos,labels,styles): leg.AddEntry(plot,label,style)
    leg.Draw()
    ## assign it to a global variable so it's not deleted
    global legend_
    legend_ = leg
    return leg


def getFriend(inputfile):
    friendfile = "{dirn}/{base}_Friend.root".format(dirn=os.path.dirname(inputfile),
                                                    base=os.path.basename(inputfile).split('.root')[0])
    return friendfile

def response_vsrun(inputfile,paramsfile,Z,panda=None):
    GBR = GBRLikelihoodTrainer(paramsfile)
    eps=3
    zSel = {'sc_truez': (Z-eps,Z+eps)}
    X,y = GBR.get_dataset(inputfile,getFriend(inputfile),loadPanda=panda,addCuts=zSel)
    filename = "gbrLikelihood_%s_mse.sav" % GBR.target.replace('sc_','')
    model = joblib.load(filename)
    y_pred = model.predict(X)
    y_raw = X[:,GBR.rawyindex]

    if GBR.target=='sc_trueint':
        YoYt = np.divide(y_pred,y)
        YrawoYt = np.divide(y_raw,y)
        xmin,xmax=0.,2.
        #print ("E/Et: ",YoYt)
    elif GBR.target=='sc_truez':
        YoYt = np.subtract(y_pred,y)
        YrawoYt = np.ones(len(y))*1e6
        xmin,xmax=-25.,25.
        #print ("Z residuals: ",YoYt," (cm)")
    else:
        RuntimeError("Target %s not foreseen." % GBR.target)


    h = ROOT.TH1F('h','',800,xmin,xmax)
    histos = {}

    histos['uncorr'] = h.Clone('h_uncorr')
    histos['regr'] = h.Clone('h_regr')

    fill_hist(histos['uncorr'],YrawoYt)
    fill_hist(histos['regr'],YoYt)

    return (histos['uncorr'],histos['regr'])
    
    
def response2D(inputfile,paramsfile,panda=None):

    friendfile = "{dirn}/{base}_Friend.root".format(dirn=os.path.dirname(inputfile),
                                                    base=os.path.basename(inputfile).split('.root')[0])
    GBR = GBRLikelihoodTrainer(paramsfile)
    X,y = GBR.get_dataset(inputfile,getFriend(inputfile),loadPanda=panda)

    xy_indices = []
    vars = GBR.variables()
    for i,v in enumerate(vars):
        if v=='sc_xmean': xy_indices.append(i)
        if v=='sc_ymean': xy_indices.append(i)

    xy = X[:,xy_indices[0]:xy_indices[1]+1]

    
    filename = "gbrLikelihood_%s_mse.sav" % GBR.target.replace('sc_','')
    model = joblib.load(filename)
    y_pred = model.predict(X)
    y_raw = X[:,GBR.rawyindex]

    if GBR.target=='sc_trueint':
        YoYt = np.divide(y_pred,y)
        YrawoYt = np.divide(y_raw,y)
        zmin,zmax=0,2
    elif GBR.target=='sc_truez':
        YoYt = np.subtract(y_pred,y)
        YrawoYt = np.ones(len(y))*1e6
        zmin,zmax=-25.,25.
    else:
        RuntimeError("Target %s not foreseen." % GBR.target)

    energy_2D = ROOT.TH3D('energy_2D','',144,0,2304,144,0,2304,200,zmin,zmax)
    energy_1D = ROOT.TH2D('energy_1D','',20,0,2304/math.sqrt(2.),70,zmin,zmax)

    energy_2Ds = {}
    energy_1Ds = {}

    energy_2Ds['uncorr'] = energy_2D.Clone('energy_2D_uncorr')
    energy_2Ds['regr'] = energy_2D.Clone('energy_2D_regr')

    energy_1Ds['uncorr'] = energy_1D.Clone('energy_1D_uncorr')
    energy_1Ds['regr'] = energy_1D.Clone('energy_1D_regr')

    hresp = ROOT.TH1F('hresp','',50,zmin,zmax)
    center = 2304/2.
    
    for i in range(len(y)):
        #print("ev {iev} has x,y=({x},{y}) and y = {z}".format(iev=i,x=xy[i][0],y=xy[i][1],z=y[i]))
        energy_2Ds['uncorr'].Fill(xy[i][1],xy[i][0],YrawoYt[i])
        energy_2Ds['regr'].Fill(xy[i][1],xy[i][0],YoYt[i])
        energy_1Ds['uncorr'].Fill(math.hypot(xy[i][1]-center,xy[i][0]-center),YrawoYt[i])
        energy_1Ds['regr'].Fill(math.hypot(xy[i][1]-center,xy[i][0]-center),YoYt[i])

    energy_2D_mode = ROOT.TH2D('energy_2D_mode','',144,0,2304,144,0,2304)
    energy_1D_mode = ROOT.TH1D('energy_1D_mode','',20,200,2100/math.sqrt(2.))
    energy_2D_mode.GetXaxis().SetTitle("ix")
    energy_2D_mode.GetYaxis().SetTitle("iy")
    energy_2D_mode.GetZaxis().SetTitle("mode")
    if GBR.target=='sc_trueint':
        energy_2D_mode.GetZaxis().SetRangeUser(0.2,1.4)
    else:
        energy_2D_mode.GetZaxis().SetRangeUser(-10,10) 
    
    energy_2D_modes = {}
    energy_2D_modes['uncorr'] = energy_2D_mode.Clone('energy_2D_mode_uncorr')
    energy_2D_modes['regr'] = energy_2D_mode.Clone('energy_2D_mode_regr')

    ROOT.gStyle.SetPaintTextFormat("1.2f");

    c = getCanvas('c')
    for k,h3D in energy_2Ds.items():
        for ix in range(1,h3D.GetNbinsX()+1):
            for iy in range(1,h3D.GetNbinsY()+1):
                # hz = ROOT.TH1F('hz','',200,0,2)
                # for iz in range(1,h3D.GetNbinsZ()+1):
                #     hz.SetBinContent(iz,h3D.GetBinContent(ix,iy,iz))
                # zmean = hz.GetMean()
                # print ("zmean = ",zmean)
                # energy_2D_modes[k].SetBinContent(ix,iy,zmean)
                                                 
                # h_integral = sum([h3D.GetBinContent(ix,iy,iz) for iz in range(1,h3D.GetNbinsZ()+1)])
                # h_runint = 0; zbinmedian = -1
                # for iz in range(1,h3D.GetNbinsZ()+1):
                #     h_runint = h_runint + h3D.GetBinContent(ix,iy,iz)
                #     if h_runint > 0.5 * h_integral:
                #         zbinmedian = iz
                #         break

                maxZ = -1; zbinmax = -1
                for iz in range(1,h3D.GetNbinsZ()+1):
                    if h3D.GetBinContent(ix,iy,iz) > maxZ:
                        maxZ = h3D.GetBinContent(ix,iy,iz)
                        zbinmax = iz
                energy_2D_modes[k].SetBinContent(ix,iy,h3D.GetZaxis().GetBinCenter(zbinmax))
        energy_2D_modes[k].Draw("colz") # text45")
        c.SaveAs("{target}_2D_{name}.png".format(target=GBR.target.replace('sc_',''),name=k))

    energy_2D_modes['ratio'] = energy_2D_modes['regr'].Clone('energy_2D_mode_ratio')
    energy_2D_modes['ratio'].Divide(energy_2D_modes['uncorr'])
    energy_2D_modes['ratio'].Draw("colz") # text45")
    c.SaveAs("%s_2D_ratio.png" % GBR.target.replace('sc_',''))
    
    
    # energy_1Ds['uncorr'].SetMarkerColor(ROOT.kBlack)
    # energy_1Ds['regr'].SetMarkerColor(ROOT.kRed)
    # energy_1Ds['uncorr'].Draw("pe1")
    # energy_1Ds['regr'].Draw("pe1 same")
    # c.SaveAs("energy_1D.png")

def makeResponseHistos(treeFile,params_txt,outfileHistos="response_histos.root",panda=None):
    outfile = ROOT.TFile.Open(outfileHistos,"recreate")
    response2D(treeFile,params_txt,panda=panda)
    outfile.cd()
    for i,dist in enumerate(Zsteps):
        print("filling graph with point ",i," and Z = ",dist)
        histos = response_vsrun(args[0],params_txt,dist,panda=panda)
        for ih,h in enumerate(histos):
            h.Draw("hist")
            suff = 'uncorr' if ih==0 else 'regr'
            name = "resp_{suf}_{dist}".format(suf=suff,dist=dist)
            h.SetName(name)
            outfile.cd()
            h.Write()                        
    outfile.Close()

def fitResponseHisto(target,histo,xmin=0.3,xmax=1.3,rebin=4,marker=ROOT.kFullCircle,color=ROOT.kRed+1):

    histo.Rebin(rebin)
    work = ROOT.RooWorkspace()
    if target=='sc_trueint':
        work.factory('Cruijff::cb(x[{xmin},{xmax}],mean[0.3,1.4],sigma[0.05,0.01,0.50],sigma,alphaL[0.1,0.01,10],alphaR[0.1,0.01,10])'.format(xmin=xmin,xmax=xmax))
    elif target=='sc_truez':
        work.factory('Cruijff::cb(x[{xmin},{xmax}],mean[-10,10],sigma[5,1,15],sigma,alphaL[0.1,0.01,10],alphaR[0.1,0.01,10])'.format(xmin=xmin,xmax=xmax))
    else:
        RuntimeError("Target %s not foreseen" % target)
    work.Print()
    
    x = work.var('x')

    histoname = histo.GetName()
    rooData = ROOT.RooDataHist(histoname,histoname,ROOT.RooArgList(work.var("x")),histo)
    getattr(work,'import')(rooData)

    frame = x.frame()
    frame.SetTitle('')

    # fit histogram
    rooData.plotOn(frame,ROOT.RooFit.MarkerStyle(marker),ROOT.RooFit.MarkerSize(2))
    pdf = work.pdf('cb')
    pdf.fitTo(rooData,ROOT.RooFit.Save())#,ROOT.RooFit.PrintLevel(-1))
    pdf.plotOn(frame,ROOT.RooFit.LineColor(color))
    rooData.plotOn(frame,ROOT.RooFit.MarkerStyle(marker),ROOT.RooFit.MarkerSize(2))

    frame.GetXaxis().SetTitleFont(42)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetTitleOffset(1.1)
    frame.GetXaxis().SetLabelFont(42)
    frame.GetXaxis().SetLabelSize(0.045)
    frame.GetXaxis().SetLabelOffset(0.007)
    frame.GetYaxis().SetTitleFont(42)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetTitleOffset(1.3)
    frame.GetYaxis().SetLabelFont(42)
    frame.GetYaxis().SetLabelSize(0.045)
    frame.GetYaxis().SetLabelOffset(0.007)
    frame.GetYaxis().SetTitle("Events")
    if target=='sc_trueint':
        frame.GetXaxis().SetTitle("E / E_{true}")
    else:
        frame.GetXaxis().SetTitle("Z - Z_{true} (cm)")        
    frame.GetXaxis().SetNdivisions(510)

    m = work.var('mean').getVal()
    s = work.var('sigma').getVal()

    em = work.var('mean').getError()
    es = work.var('sigma').getError()

    return (frame,m,s,em,es)

def fitAllResponses(inputfile="response_histos.root",target='sc_trueint'):

    if target=='sc_trueint':
        limits = {5 : (0.1,1.6),
                  15  : (0.3,1.5),
                  25  : (0.3,1.5),
                  36  : (0.3,1.5),
                  48  : (0.4,1.5)}
    elif target=='sc_truez':
        limits = {5 : (-25,25),
                  15  : (-25,25),
                  25  : (-25,25),
                  36  : (-25,25),
                  48  : (-25,25)}

    dummy_uncorr = ROOT.TH1F("dummy_uncorr","",1,0,1)
    dummy_regr = ROOT.TH1F("dummy_regr","",1,0,1)

    dummy_uncorr.SetMarkerStyle(ROOT.kOpenSquare)
    dummy_uncorr.SetLineColor(ROOT.kRed+2)
    dummy_uncorr.SetMarkerSize(2)
    
    dummy_regr.SetMarkerStyle(ROOT.kFullCircle)
    dummy_regr.SetLineColor(ROOT.kBlue)
    dummy_regr.SetMarkerSize(2)

    results = {}
    tf = ROOT.TFile(inputfile,"read")
    for i,dist in enumerate(Zsteps):
        c = getCanvas()
        corrs = ["regr","uncorr"] if target=='sc_trueint' else ['regr']
        for corr in corrs:
            hname = "resp_{corr}_{dist}".format(corr=corr,dist=dist)
            histo = tf.Get(hname)
            xmin = limits[dist][0]; xmax=limits[dist][1]
            if corr=="uncorr":
                marker = ROOT.kOpenSquare
                linecol = ROOT.kRed+2
            else:
                marker = ROOT.kFullCircle
                linecol = ROOT.kBlue
            (frame,mean,sigma,meanerr,sigmaerr)=fitResponseHisto(target,histo,xmin,xmax,8,marker,linecol)
            frame.Draw("same")
            results[(corr,dist)] = (mean,sigma,meanerr,sigmaerr)
            
        if target=='sc_trueint':
            responses = [dummy_uncorr,dummy_regr]
            titles = ['E_{raw}','E_{regr}']
            styles = ['pl','pl']
            legend = doLegend(responses,titles,styles,corner="TL")
        else:
            responses = [dummy_regr]
            titles = ['Z_{regr}']
            styles = ['pl']
            legend = doLegend(responses,titles,styles,corner="TL")

        doTinyCmsPrelim(hasExpo = False,textSize=0.04, options=None,doWide=False)
        for ext in ['pdf','png','root']:
            c.SaveAs('respcomp_{target}_{dist}_fit.{ext}'.format(target=target,dist=dist,ext=ext))
    return results

def graphVsR(typeInput='histo',histoFile="response_histos.root",target='sc_trueint'):

    resp_vs_z_uncorr = ROOT.TGraphErrors(len(Zsteps))
    reso_vs_z_uncorr = ROOT.TGraph(len(Zsteps))
    resp_vs_z_regr = ROOT.TGraphErrors(len(Zsteps))
    reso_vs_z_regr = ROOT.TGraph(len(Zsteps))

    reso_vs_z_uncorr_fullrms = ROOT.TGraph(len(Zsteps))
    reso_vs_z_regr_fullrms = ROOT.TGraph(len(Zsteps))

    resp_vs_z_regr.GetXaxis().SetTitle("z (cm)")
    reso_vs_z_regr.GetXaxis().SetTitle("z (cm)")    
    resp_vs_z_regr.SetTitle("")
    reso_vs_z_regr.SetTitle("")
    if target=='sc_trueint':    
        resp_vs_z_regr.GetYaxis().SetTitle("E / E_{true}")
        reso_vs_z_regr.GetYaxis().SetTitle("E / E_{true} resolution")
    else:
        resp_vs_z_regr.GetYaxis().SetTitle("Z - Z_{true} (cm)")
        reso_vs_z_regr.GetYaxis().SetTitle("Z resolution (cm)")

    
    if typeInput=="histo":
        tf = ROOT.TFile(histoFile,"read")
        for i,dist in enumerate(Zsteps):
            huncorr = tf.Get("resp_uncorr_{dist}".format(dist=dist))
            hregr   = tf.Get("resp_regr_{dist}".format(dist=dist))
            resp_vs_z_uncorr.SetPoint(i,dist,huncorr.GetMean())
            resp_vs_z_uncorr.SetPointError(i,0,huncorr.GetMeanError())

            resp_vs_z_regr.SetPoint(i,dist,hregr.GetMean())
            resp_vs_z_regr.SetPointError(i,0,hregr.GetMeanError())
        
            reso_vs_z_uncorr.SetPoint(i,dist,huncorr.GetRMS()/huncorr.GetMean())
            reso_vs_z_regr.SetPoint(i,dist,hregr.GetRMS()/hregr.GetMean())
        tf.Close()

    elif typeInput=="fit":
        results = fitAllResponses(target=target)
        tf = ROOT.TFile(histoFile,"read")
        for i,dist in enumerate(Zsteps):

            if target=='sc_trueint':
                mean,sigma,meanerr,sigmaerr = results[('uncorr',dist)]
                resp_vs_z_uncorr.SetPoint(i,dist,mean)
                resp_vs_z_uncorr.SetPointError(i,0,meanerr)
                reso_vs_z_uncorr.SetPoint(i,dist,sigma/mean)
                huncorr = tf.Get("resp_uncorr_{dist}".format(dist=dist))
                reso_vs_z_uncorr_fullrms.SetPoint(i,dist,huncorr.GetRMS()/huncorr.GetMean())

            mean,sigma,meanerr,sigmaerr = results[('regr',dist)]
            resp_vs_z_regr.SetPoint(i,dist,mean)
            resp_vs_z_regr.SetPointError(i,0,meanerr)
            hregr   = tf.Get("resp_regr_{dist}".format(dist=dist))
            if target=='sc_trueint':
                reso_vs_z_regr.SetPoint(i,dist,sigma/mean)
                reso_vs_z_regr_fullrms.SetPoint(i,dist,hregr.GetRMS()/hregr.GetMean())
            else:
                reso_vs_z_regr.SetPoint(i,dist,sigma)                
                reso_vs_z_regr_fullrms.SetPoint(i,dist,hregr.GetRMS())

        tf.Close()

    c = getCanvas('c')

    # response
    resp_vs_z_uncorr.SetMarkerStyle(ROOT.kFullCircle)
    resp_vs_z_regr.SetMarkerStyle(ROOT.kFullCircle)
    resp_vs_z_uncorr.SetMarkerSize(2)
    resp_vs_z_regr.SetMarkerSize(2)
    resp_vs_z_uncorr.SetMarkerColor(ROOT.kRed)
    resp_vs_z_uncorr.SetLineColor(ROOT.kRed)
    resp_vs_z_regr.SetMarkerColor(ROOT.kBlack)
    resp_vs_z_regr.Draw("Ape")
    if target=='sc_trueint':
        resp_vs_z_uncorr.Draw("pe")

    if target=='sc_trueint':
        responses = [resp_vs_z_uncorr,resp_vs_z_regr]
        titles = ['raw','regression']
        styles = ['pl','pl']
    else:
        responses = [resp_vs_z_regr]
        titles = ['regression']
        styles = ['pl']        
    
    if target=='sc_trueint':
        resp_vs_z_regr.GetYaxis().SetRangeUser(0.2,1.3)
    else:
        resp_vs_z_regr.GetYaxis().SetRangeUser(-10,10)
    legend = doLegend(responses,titles,styles,corner="BR")    
    for ext in ['pdf','png']:
        c.SaveAs("response_%s.%s"%(target.replace('sc_',''),ext))

    # resolution
    reso_vs_z_uncorr.SetMarkerStyle(ROOT.kFullCircle)
    reso_vs_z_regr.SetMarkerStyle(ROOT.kFullCircle)
    reso_vs_z_uncorr.SetMarkerColor(ROOT.kRed)
    reso_vs_z_regr.SetMarkerColor(ROOT.kBlack)
    reso_vs_z_uncorr.SetMarkerSize(2)
    reso_vs_z_regr.SetMarkerSize(2)
    reso_vs_z_regr.Draw("Ape")
    if target=='sc_trueint':
        reso_vs_z_uncorr.Draw("pe")
    legCols = 1
    if typeInput=='fit':
        reso_vs_z_uncorr_fullrms.SetMarkerStyle(ROOT.kOpenSquare)
        reso_vs_z_regr_fullrms.SetMarkerStyle(ROOT.kOpenSquare)
        reso_vs_z_uncorr_fullrms.SetMarkerColor(ROOT.kRed)
        reso_vs_z_regr_fullrms.SetMarkerColor(ROOT.kBlack)
        reso_vs_z_uncorr_fullrms.SetMarkerSize(2)
        reso_vs_z_regr_fullrms.SetMarkerSize(2)
        if target=='sc_trueint':
            responses = [resp_vs_z_uncorr,resp_vs_z_regr,reso_vs_z_uncorr_fullrms,reso_vs_z_regr_fullrms]
            titles = ['uncorrected (#sigma_{G})','regression (#sigma_{G})','uncorrected (rms)','regression (rms)']
            styles = ['pl','pl','pl','pl']
            reso_vs_z_uncorr_fullrms.Draw("pe")
        else:
            responses = [resp_vs_z_regr,reso_vs_z_regr_fullrms]            
            titles = ['regression (#sigma_{G})','regression (rms)']
            styles = ['pl','pl']

        reso_vs_z_regr_fullrms.Draw("pe")
        legCols=2

    if target=='sc_trueint':
        reso_vs_z_regr.GetYaxis().SetRangeUser(0.0,0.4)
    else:
        reso_vs_z_regr.GetYaxis().SetRangeUser(0,15)
    legend = doLegend(responses,titles,styles,corner="TL",textSize=0.03,legWidth=0.8,nColumns=legCols)
    for ext in ['pdf','png']:
        c.SaveAs("resolution_%s.%s"%(target.replace('sc_',''),ext))

            
if __name__ == '__main__':

    from optparse import OptionParser
    parser = OptionParser(usage='%prog input.root params_gbrtrain.txt [opts] ')
    parser.add_option(        '--loadPanda', dest='loadPanda', default=None, type='string', help='file with regression data as panda dataframe to be loaded instead of the full ROOT files')
    parser.add_option(        '--target', dest='target', default="sc_trueint", type='string', help='target variable of the regression (can be sc_trueint for energy or sc_truez for z)')
    (options, args) = parser.parse_args()

    makeResponseHistos(args[0],args[1],panda=options.loadPanda)
    graphVsR("fit",target=options.target)
