import ROOT
import uproot
import numpy as np
import pandas as pd
import glob
import re
import tqdm
from multiprocessing import Pool
import argparse
import dask
import dask.array as da
from dask import delayed

def get_numbers_from_filename(filename):
    return re.search(r'\d+', filename).group(0)
def grapherr(x,y,ex,ey,x_string, y_string,name=None, color=9, markerstyle=33, markersize=2,write=True):
        plot = ROOT.TGraphErrors(len(x),  np.array(x  ,dtype="d")  ,   np.array(y  ,dtype="d") , np.array(   ex   ,dtype="d"),np.array( ey   ,dtype="d"))
        if name is None: plot.SetNameTitle(y_string+" vs "+x_string,y_string+" vs "+x_string)
        else: plot.SetNameTitle(name, name)
        plot.GetXaxis().SetTitle(x_string)
        plot.GetYaxis().SetTitle(y_string)
        plot.SetMarkerColor(color)#blue
        plot.SetMarkerStyle(markerstyle)
        plot.SetMarkerSize(markersize)
        if write==True: plot.Write()
        return plot
def fill_h(histo_name, array):
    for x in range (len(array)):
        histo_name.Fill((np.array(array[x] ,dtype="d")))
def hist(list, x_name, channels=1000, linecolor=4, linewidth=4,write=True,startZero=False,BinSpacing=None):
    array=np.array(list ,dtype="d")
    if BinSpacing is None: ch=channels
    else: ch = int((np.max(array) - np.min(array)) / BinSpacing)
    if startZero==False: hist=ROOT.TH1D(x_name,x_name,ch,0.99*np.min(array),1.01*np.max(array))
    else: hist=ROOT.TH1D(x_name,x_name,ch,0,1.01*np.max(array))
    fill_h(hist,array)
    hist.SetLineColor(linecolor)
    hist.SetLineWidth(linewidth)
    hist.GetXaxis().SetTitle(x_name)
    hist.GetYaxis().SetTitle("Entries")
    if write==True: hist.Write()
    #hist.SetStats(False)
    hist.GetYaxis().SetMaxDigits(3);
    hist.GetXaxis().SetMaxDigits(3);
    return hist
def process_file(file):
    run_number = get_numbers_from_filename(file)
    root_file = uproot.open(file)

    mean_h2 = root_file["pedmap"]
    rms_h2 = root_file["pedmapsigma"]

    # Create Dask arrays from Numpy arrays
    flat_means = da.from_array(mean_h2.to_numpy()[0].flatten(), chunks=1000)  # Adjust the chunk size as needed
    flat_rms = da.from_array(rms_h2.to_numpy()[0].flatten(), chunks=1000)

    # Perform calculations on Dask arrays
    noisyPix110 = da.count_nonzero(flat_means >= 110).compute()
    noisyPix103 = da.count_nonzero(flat_means >= 103).compute()
    noisyPix97 = da.count_nonzero(flat_means <= 97).compute()
    mean_value = flat_means.mean().compute()
    integral_value = flat_means.sum().compute()

    rmsnoisyPix5 = da.count_nonzero(flat_rms >= 5).compute()
    rmsnoisyPix10 = da.count_nonzero(flat_rms >= 10).compute()
    rmsnoisyPix15 = da.count_nonzero(flat_rms >= 15).compute()
    rms_value = flat_rms.mean().compute()

    return run_number, noisyPix110, noisyPix103, noisyPix97, rmsnoisyPix5, rmsnoisyPix10, rmsnoisyPix15, mean_value, rms_value, integral_value

#Parser
parser = argparse.ArgumentParser(description='Analyze pedestal', epilog='Version: 1.0')
parser.add_argument('-j','--cores',help='select # cores', action='store', type=int, default=1)
parser.add_argument('-t','--test',help='select test mode', action='store', type=bool, default=False)
args = parser.parse_args()
# Number of concurrent processes
num_processes = args.cores
#check all pedestals
if args.test==True: files=glob.glob("ped_test/*.root")
else: files=glob.glob("pedestals/*.root")
run_num,noisyPix110,noisyPix103,noisyPix97,rmsnoisyPix5,rmsnoisyPix10,rmsnoisyPix15,mean,rms,integral=np.empty(len(files)),np.empty(len(files)),np.empty(len(files)),np.empty(len(files)),np.empty(len(files)),np.empty(len(files)),np.empty(len(files)),np.empty(len(files)),np.empty(len(files)),np.empty(len(files))
files.reverse()
#create Pool
with Pool(num_processes) as pool:
    results = list(tqdm.tqdm(pool.imap(process_file, files), total=len(files)))
## Unpack the results
for i, (run_num_result, noisyPix110_result, noisyPix103_result, noisyPix97_result, rmsnoisyPix5_result, rmsnoisyPix10_result, rmsnoisyPix15_result, mean_result, rms_result, integral_result) in enumerate(results):
    run_num[i] = run_num_result
    noisyPix110[i] = noisyPix110_result
    noisyPix103[i] = noisyPix103_result
    noisyPix97[i] = noisyPix97_result
    mean[i] = mean_result
    integral[i] = integral_result
    rmsnoisyPix5[i] = rmsnoisyPix5_result
    rmsnoisyPix10[i] = rmsnoisyPix10_result
    rmsnoisyPix15[i] = rmsnoisyPix15_result
    rms[i] = rms_result

    # for example
    if i == (len(files) - 1):
        mean_hist = uproot.open(files[i])["pedmean"]
        rms_hist = uproot.open(files[i])["pedrms"]
#output file root file creation
main=ROOT.TFile("anal_ped.root","RECREATE")
"""
#COMPARISON REAL AND COMPUTED RMS ADN MEAN
#not good in multicore
hist(flat_means,"pedmean Example Run "+str(run_num[-1]),BinSpacing=0.01)
hist(flat_rms,"pedrms Example Run "+str(run_num[-1]),startZero=True,BinSpacing=0.01)
"""
# Create a ROOT TH1F histogram rms from pedfile
bin_edges = rms_hist.to_numpy()[1]
bin_entries = rms_hist.to_numpy()[0]
hist_rms_pedfile = ROOT.TH1F("rms from Pedfile Run "+str(run_num[-1]), "rms from Pedfile Run "+str(run_num[-1]), len(bin_entries), bin_edges)
for i,entry in enumerate(bin_entries):
    hist_rms_pedfile.SetBinContent(i+1, entry)
hist_rms_pedfile.Write()
# Create a ROOT TH1F histogram mean from pedfile
bin_edges = mean_hist.to_numpy()[1]
bin_entries = mean_hist.to_numpy()[0]
hist_mean_pedfile = ROOT.TH1F("mean from Pedfile Run "+str(run_num[-1]), "mean from Pedfile Run "+str(run_num[-1]), len(bin_entries), bin_edges)
for i,entry in enumerate(bin_entries):
    hist_mean_pedfile.SetBinContent(i+1, entry)
hist_mean_pedfile.Write()

#CHECK VARIABLES
def canvas_history(x,y,x_string,y_string,name):
    graph=grapherr(x,y,1E-20*np.ones(len(run_num)),1E-20*np.ones(len(run_num)),x_string,y_string,name)
    can1=ROOT.TCanvas(name, name, 1200   ,1200)
    can1.SetFillColor(0)
    can1.SetBorderMode(0)
    can1.SetBorderSize(2)
    can1.SetLeftMargin(0.18)
    can1.SetRightMargin(0.02)
    can1.SetTopMargin(0.1)
    can1.SetBottomMargin(0.1)
    can1.SetFrameBorderMode(0)
    can1.SetFrameBorderMode(0)
    can1.SetFixedAspectRatio()
    #can1.cd()
    graph.Draw("AP")
    #graph.GetXaxis().SetDecimals(1)
    #graph.GetXaxis().SetMaxDigits(2)

    can1.Update()
    ymax=ROOT.gPad.GetUymax()
    ymin=ROOT.gPad.GetUymin()
    xmax=ROOT.gPad.GetUxmax()
    xmin=ROOT.gPad.GetUxmin()

    #lines 23802,24328,AmBe Campaign with also GEM Off
    line1=ROOT.TLine(23.802,ymin,23.802,ymax)
    line1.SetLineColor(40)
    line1.SetLineWidth(2)
    line1.Draw("SAME")
    line2=ROOT.TLine(25.427,ymin,25.427,ymax)
    line2.SetLineColor(40)
    line2.SetLineWidth(2)
    line2.Draw("SAME")
    pave=ROOT.TPaveText(23.802,ymax-0.02*(ymax-ymin),25.427,ymax-0.02*(ymax-ymin))
    pave.SetBorderSize(2)
    pave.SetFillColor(0)
    pave.SetTextSize(0.02)
    pave.SetTextColor(40)
    pave.AddText("AmBe Campaign")
    pave.Draw("SAME")

    #lines 27858,30831,Ba Campaign with also collimator off
    line3=ROOT.TLine(27.858,ymin,27.858,ymax)
    line3.SetLineColor(41)
    line3.SetLineWidth(2)
    line3.Draw("SAME")
    line4=ROOT.TLine(30.831,ymin,30.831,ymax)
    line4.SetLineColor(41)
    line4.SetLineWidth(2)
    line4.Draw("SAME")
    pave1=ROOT.TPaveText(27.858,ymax-0.02*(ymax-ymin),30.831,ymax-0.02*(ymax-ymin))
    pave1.SetBorderSize(2)
    pave1.SetFillColor(0)
    pave1.SetTextSize(0.02)
    pave1.SetTextColor(41)
    pave1.AddText("Ba Campaign")
    pave1.Draw("SAME")

    #lines 30845,36663,Ba Campaign with also collimator off
    line5=ROOT.TLine(30.845,ymin,30.845,ymax)
    line5.SetLineColor(42)
    line5.SetLineWidth(2)
    line5.Draw("SAME")
    line6=ROOT.TLine(36.663,ymin,36.663,ymax)
    line6.SetLineColor(42)
    line6.SetLineWidth(2)
    line6.Draw("SAME")
    pave2=ROOT.TPaveText(30.845,ymax-0.02*(ymax-ymin),36.663,ymax-0.02*(ymax-ymin))
    pave2.SetBorderSize(2)
    pave2.SetFillColor(0)
    pave2.SetTextSize(0.02)
    pave2.SetTextColor(42)
    pave2.AddText("Eu Campaign")
    pave2.Draw("SAME")

    #lines 37955,38292,Am Campaign
    line7=ROOT.TLine(37.955,ymin,37.955,ymax)
    line7.SetLineColor(43)
    line7.SetLineWidth(2)
    line7.Draw("SAME")
    line8=ROOT.TLine(38.292,ymin,38.292,ymax)
    line8.SetLineColor(43)
    line8.SetLineWidth(2)
    line8.Draw("SAME")
    pave3=ROOT.TPaveText(37.955,ymax-0.02*(ymax-ymin),38.292,ymax-0.02*(ymax-ymin))
    pave3.SetBorderSize(2)
    pave3.SetFillColor(0)
    pave3.SetTextSize(0.02)
    pave3.SetTextColor(43)
    pave3.AddText("Am Campaign")
    pave3.Draw("SAME")

    can1.Write()
    can1.SaveAs("output_plots/"+name+".png")

def canvas_norm(x,y,x_string,y_string,name):
    graph=grapherr(x,y,1E-20*np.ones(len(run_num)),1E-20*np.ones(len(run_num)),x_string,y_string,name)
    can1=ROOT.TCanvas(name, name, 1200   ,1200)
    can1.SetFillColor(0)
    can1.SetBorderMode(0)
    can1.SetBorderSize(2)
    can1.SetLeftMargin(0.18)
    can1.SetRightMargin(0.02)
    can1.SetTopMargin(0.1)
    can1.SetBottomMargin(0.1)
    can1.SetFrameBorderMode(0)
    can1.SetFrameBorderMode(0)
    can1.SetFixedAspectRatio()
    #can1.cd()
    graph.Draw("AP")
    #graph.GetXaxis().SetDecimals(1)
    #graph.GetXaxis().SetMaxDigits(2)

    can1.Write()
    can1.SaveAs("output_plots/"+name+".png")


#ONLY GRAPHERR
canvas_norm(run_num/1000,mean,"Pedestal Run Number (#times 10^{3})","Mean from TH2","Mean from TH2")
canvas_norm(run_num/1000,rms,"Pedestal Run Number (#times 10^{3})","RMS from TH2","RMS from TH2")
canvas_norm(run_num/1000,integral,"Pedestal Run Number (#times 10^{3})","Integral from TH2","Integral from TH2")
canvas_norm(run_num/1000,mean/np.max(mean),"Pedestal Run Number (#times 10^{3})","Normalized Mean from TH2","Normalized Mean from TH2")
canvas_norm(run_num/1000,rms/np.max(rms),"Pedestal Run Number (#times 10^{3})","Normalized RMS from TH2","Normalized RMS from TH2")
canvas_norm(run_num/1000,integral/np.max(integral),"Pedestal Run Number (#times 10^{3})","Normalized Integral from TH2","Normalized Integral from TH2")



"""
canvas_history(run_num/1000,noisyPix110,"Pedestal Run Number (#times 10^{3})","# Pixels mean>110","# Pixels mean>110")
canvas_history(run_num/1000,noisyPix103,"Pedestal Run Number (#times 10^{3})","# Pixels mean>103","# Pixels mean>103")
canvas_history(run_num/1000,noisyPix97,"Pedestal Run Number (#times 10^{3})","# Pixels mean<97","# Pixels mean<97")
canvas_history(run_num/1000,rmsnoisyPix5,"Pedestal Run Number (#times 10^{3})","# Pixels rms>5","# Pixel rms>5")
canvas_history(run_num/1000,rmsnoisyPix10,"Pedestal Run Number (#times 10^{3})","# Pixels rms>10","# Pixel rms>10")
canvas_history(run_num/1000,rmsnoisyPix15,"Pedestal Run Number (#times 10^{3})","# Pixels rms>15","# Pixel rms>15")
canvas_history(run_num/1000,mean,"Pedestal Run Number (#times 10^{3})","Mean from TH2","Mean from TH2")
canvas_history(run_num/1000,rms,"Pedestal Run Number (#times 10^{3})","RMS from TH2","RMS from TH2")
canvas_history(run_num/1000,integral,"Pedestal Run Number (#times 10^{3})","Integral from TH2","Integral from TH2")

#normalized
canvas_history(run_num/1000,noisyPix110/np.max(noisyPix110),"Pedestal Run Number (#times 10^{3})","Norm # Pixels mean>110","Norm # Pixels mean>110")
canvas_history(run_num/1000,noisyPix103/np.max(noisyPix103),"Pedestal Run Number (#times 10^{3})","Norm # Pixels mean>103","Norm # Pixels mean>103")
canvas_history(run_num/1000,noisyPix97/np.max(noisyPix97),"Pedestal Run Number (#times 10^{3})","Norm # Pixels mean<97","Norm # Pixels mean<97")
canvas_history(run_num/1000,rmsnoisyPix5/np.max(rmsnoisyPix5),"Pedestal Run Number (#times 10^{3})","Norm # Pixels rms>5","Norm # Pixel rms>5")
canvas_history(run_num/1000,rmsnoisyPix10/np.max(rmsnoisyPix10),"Pedestal Run Number (#times 10^{3})","Norm # Pixels rms>10","Norm # Pixel rms>10")
canvas_history(run_num/1000,rmsnoisyPix15/np.max(rmsnoisyPix15),"Pedestal Run Number (#times 10^{3})","Norm # Pixels rms>15","Norm # Pixel rms>15")
canvas_history(run_num/1000,mean/np.max(mean),"Pedestal Run Number (#times 10^{3})","Norm Mean from TH2","Norm Mean from TH2")
canvas_history(run_num/1000,rms/np.max(rms),"Pedestal Run Number (#times 10^{3})","Norm RMS from TH2","Norm RMS from TH2")
canvas_history(run_num/1000,integral/np.max(integral),"Pedestal Run Number (#times 10^{3})","Norm Integral from TH2","Norm Integral from TH2")

"""