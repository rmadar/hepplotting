import ROOT
import numpy as np
import pandas as pd
import math
import sys
import os

ROOT.gROOT.LoadMacro(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'AtlasStyle.C'))
ROOT.SetAtlasStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)


def ATLASLabel(x, y, text, withRatio=True, rsize=None):
    txt = ROOT.TLatex()
    ROOT.SetOwnership(txt, False)
    txt.SetNDC()
    delx, size = 0.132, 0.052
    if rsize:
        delx *= rsize   # To account for not 800x700 canvas
    if withRatio:
        size = 0.068
    txt.SetTextSize(size)
    txt.SetTextFont(72)
    txt.DrawLatex(x, y, 'ATLAS')
    txt.SetTextFont(42)
    txt.DrawLatex(x+delx, y, text)


def stampText(text, x, y, size):
    txt = ROOT.TLatex()
    ROOT.SetOwnership(txt, False)
    txt.SetNDC()
    txt.SetTextFont(42)
    txt.SetTextColor(1)
    txt.SetTextSize(size)
    txt.DrawLatex(x, y, text)


def sum_histograms(hBkg, name='tot'):
    '''
     Histogram summer
     ================

    - Args:
    . "hBkg" [list of TH1] is a list of histograms to be summed up

    - Return:
    . TH1 being the summed histogram
    '''
    for i, h in enumerate(hBkg):
        if i == 0:
            hTot = h.Clone(h.GetName()+'_'+name)
        else:
            hTot.Add(h)
    return hTot.Clone(name)


def add_flat_syst(h, s=0, name='wsyst'):
    ''' Adding a flat systematic to a given histogram
    - Args:
    . h [TH1] is the histogram
    . s [float] is the relative systematics (ie 15% would be s=0.15)
    - Return:
    . TH1 with the inflated error (stat ++ syst)
    '''
    hres = h.Clone(h.GetName()+'_'+name)
    if s == 0:
        return hres
    else:
        for i in range(0, h.GetNbinsX()+1):
            val, err = h.GetBinContent(i), h.GetBinError(i)
            new_err = math.sqrt(err**2 + (val*s)**2)
            hres.SetBinContent(i, val)
            hres.SetBinError(i, new_err)
    return hres


def scale_xaxis(h, scale, addOverflow=False):
    '''
    Re-scale the x-axis of h [TH1] by a factor scale [double] and possibly
    add the overflow bin in the last one
    '''
    nbins = h.GetNbinsX()
    xmin = h.GetBinLowEdge(1)
    xmax = h.GetBinLowEdge(nbins)+h.GetBinWidth(nbins)
    hres = ROOT.TH1F(h.GetName()+'_goodbin', h.GetTitle(), nbins, xmin*scale, xmax*scale)
    for i in xrange(0, nbins+2):
        hres.SetBinContent(i, h.GetBinContent(i))
        hres.SetBinError(i, h.GetBinError(i))
    if addOverflow:
        Nlb, Elb = hres.GetBinContent(nbins), hres.GetBinError(nbins)
        Nof, Eof = hres.GetBinContent(nbins+1), hres.GetBinError(nbins+1)
        hres.SetBinContent(nbins, Nlb+Nof)
        hres.SetBinError(nbins, np.sqrt(Elb**2+Eof**2))
        hres.SetBinContent(nbins+1, 0), hres.SetBinError(nbins+1, 0)
    return hres


def remove_0entry_data(hdata, th=0.5):
    '''
    Set bin content and error to a not visible values for data histogram \'hdata\' [TH1].
    The threshold can be tuned (e.g. ratio to data wher undesirable values are not only <0.5).
    '''
    for i in xrange(1, hdata.GetNbinsX()+2):
        if hdata.GetBinContent(i) < th:
            hdata.SetBinContent(i, -1e5)
            hdata.SetBinError(i, 0.0)
    return hdata


def make_nice_canvas(dictBkg, hTot, hData, plot_name, **kwargs):
    '''
    Produce a canvas with stacked histograms for background, data and ratio plots.


    Required arguments
    ==================
    . dictBkg [dictionnary {bkgName:[TH1,color,legName]}] containing all background histograms
    . hTot [TH1] is the histogram of the total data with its uncertainty
    . hData [TH1] is the histogram of data
    . plot_name [string] is the name of the final plot (plot_name.pdf)


    Key-word arguments
    ==================
    . plotdir [string] is a directory where the plots will be stored (default: 'plots')

    AXIS properties
    ---------------
    . xtitle [string] is x-axis title
    . ytitle [string] is y-axis title
    . dictSig [dict {sigName:[TH1,color,norm,legName]}] is dictionnary with name [string], histo [TH1],
     color [int], norm [float] and legName [string] of several signals
    . is_logy [boolean] to plot in log scale or not
    . bin_label [list of string] to name bins (e.g plots with one region yield per bin)
    . xlabel_size [float] size of the x-axis bin labels
    . xlabel_offset [float] offset of the x-axis bin labels
    . xticksInt [bool] keep only integer values for x-axis ticks
    . xmin [float] lower x-axis value
    . xmax [float] lower x-axis value
    . ymin [float] lower y-axis value
    . ymax [float] higher y-axis value
    . r_ymin [float] lower y-axis value on the ratio plot
    . r_ymax [float] higher y-axis value on the ratio plot

    LEGEND properties
    -----------------
    . leg_pos [list of float] specify the legend position via bottom left (x1,y1)
     and top right (x2,y2) using [x1,y1,x2,y2]
    . unc_leg [string] to tune the name of uncertainty (eg. stat-only)
    . leg_ncols [int] number of columns used for the legend
    . leg_put_nevts [bool] to print events yields in the legend (default: False)
    . leg_textsize [float] size of the legend text (~0.030 to ~0.045)

    HISTO properties
    ----------------
    . m_size [float] is the marker size for data
    . error_fill [int] is the filling style for the uncertainty band
    . error_alpha [float] is the transparency for the uncertainty band (in [0,1])
    . histo_border [int] is the border size of background histograms in the stacks

    LABELS properties
    -----------------
    . plot_labels [list of string] given the labels printed below ATLAS and Lumi
    . atlas_label [string] is \'Internal\' by default but can be \'ATLAS\', \'Preliminary\', \'Simulation\'
    . lumi [float] is the integrated luminosity (default: 1/fb)

    CANVAS properties
    -----------------
    . canvas [TCanvas] on which the plot will be made
    . can_ratio [float] specify the canvas size such as width=900/ratio and height=800
    . can_scale [float] scale the whole canvas without changin its ratio
    . plot_ratio [boolean] to plot or not the ratio panel
    . ratio_type [string] to choose what to plot in the bottom plot (\'ratio\' [default], \'SoverB\', \'signif\')
    '''

    plotdir, dictSig, sig_line_style, xtitle_arg, ytitle_arg = 'plots', None, 1, None, None
    ymin_arg, ymax_arg, xmin_arg, xmax_arg, leg_pos, lumi = None, None, None, None, None, 1.0
    is_logy, bin_label, xlabel_size, xlabel_offset, xticksInt = None, None, None, None, False
    r_ymin, r_ymax, can_ratio, can_scale, m_size = None, None, None, 1.0, None
    canvas, error_fill, error_alpha, histo_border, plot_labels = None, 3356, 0.3, 0, None
    plot_ratio, atlas_label, unc_leg, ratio_type = True, 'Internal', 'Total bkg w/ unc.', 'ratio'
    leg_with_nevts, leg_ncols, leg_textsize = False, 1, None
    if 'lumi' in kwargs:
        lumi = kwargs['lumi']
    if 'dictSig' in kwargs:
        dictSig = kwargs['dictSig']
    if 'sig_line_style'in kwargs:
        sig_line_style = kwargs['sig_line_style']
    if 'plotdir' in kwargs:
        plotdir = kwargs['plotdir']
    if 'xtitle' in kwargs:
        xtitle_arg = kwargs['xtitle']
    if 'ytitle' in kwargs:
        ytitle_arg = kwargs['ytitle']
    if 'bkg_color' in kwargs:
        bkg_color = kwargs['bkg_color']
    if 'is_logy' in kwargs:
        is_logy = kwargs['is_logy']
    if 'bin_label' in kwargs:
        bin_label = kwargs['bin_label']
    if 'xlabel_size' in kwargs:
        xlabel_size = kwargs['xlabel_size']
    if 'xlabel_offset' in kwargs:
        xlabel_offset = kwargs['xlabel_offset']
    if 'xticksInt' in kwargs:
        xticksInt = kwargs['xticksInt']
    if 'xmin' in kwargs:
        xmin_arg = kwargs['xmin']
    if 'xmax' in kwargs:
        xmax_arg = kwargs['xmax']
    if 'ymin' in kwargs:
        ymin_arg = kwargs['ymin']
    if 'ymax' in kwargs:
        ymax_arg = kwargs['ymax']
    if 'r_ymin' in kwargs:
        r_ymin = kwargs['r_ymin']
    if 'r_ymax' in kwargs:
        r_ymax = kwargs['r_ymax']
    if 'canvas' in kwargs:
        canvas = kwargs['canvas']
    if 'can_ratio' in kwargs:
        can_ratio = kwargs['can_ratio']
    if 'can_scale' in kwargs:
        can_scale = kwargs['can_scale']
    if 'leg_pos' in kwargs:
        leg_pos = kwargs['leg_pos']
    if 'leg_ncols' in kwargs:
        leg_ncols = kwargs['leg_ncols']
    if 'leg_put_nevts' in kwargs:
        leg_put_nevts = kwargs['leg_put_nevts']
    if 'leg_textsize' in kwargs:
        leg_textsize = kwargs['leg_textsize']
    if 'unc_leg' in kwargs:
        unc_leg = kwargs['unc_leg']
    if 'm_size' in kwargs:
        m_size = kwargs['m_size']
    if 'plot_labels' in kwargs:
        plot_labels = kwargs['plot_labels']
    if 'atlas_label' in kwargs:
        atlas_label = kwargs['atlas_label']
    if 'error_fill' in kwargs:
        error_fill = kwargs['error_fill']
    if 'error_alpha' in kwargs:
        error_alpha = kwargs['error_alpha']
    if 'histo_border' in kwargs:
        histo_border = kwargs['histo_border']
    if 'plot_ratio' in kwargs:
        plot_ratio = kwargs['plot_ratio']
    if 'ratio_type' in kwargs:
        ratio_type = kwargs['ratio_type']

    # Get color and names for bkg histograms
    bkg_name = dictBkg.keys()
    hBkg = {n: v[0] for n, v in dictBkg.items()}
    bkg_color = {n: v[1] for n, v in dictBkg.items()}
    bkg_legname = {n: v[2] for n, v in dictBkg.items()}

    # Histo cosmetics
    if dictSig:
        for n, sig in dictSig.items():
            h, color, norm, legName = sig
            ROOT.SetOwnership(h, False)
            h.SetTitle('')
            h.SetFillColor(0)
            h.SetLineColor(color)
            h.SetLineWidth(4)
            h.SetLineStyle(sig_line_style)
            h.SetFillColor(0)
            if norm:
                if h.Integral(-999, 999)>0:
                    h.Scale(norm/h.Integral(-999, 999))
                else:
                    h.Scale(0.)
            h.SetFillColor(0)
    for b in bkg_name:
        ROOT.SetOwnership(hBkg[b], False)
        hBkg[b].SetTitle('')
        hBkg[b].SetLineWidth(histo_border)
        hBkg[b].SetMarkerSize(0)
        hBkg[b].SetFillColor(bkg_color[b])
        hBkg[b].SetLineColorAlpha(1, 0.3)
        hBkg[b].SetMarkerColor(bkg_color[b])

    ROOT.SetOwnership(hTot, False)
    hTot.SetLineWidth(0)

    # Preparing the stack
    hstack = ROOT.THStack()
    ROOT.SetOwnership(hstack, False)
    for b in bkg_name[::-1]:
        hstack.Add(hBkg[b])

    # Manage axis scaling and label involving numbers
    ROOT.SetOwnership(hData, False)
    nbins = hData.GetNbinsX()
    xmin = hData.GetBinLowEdge(1)
    xmax = hData.GetBinLowEdge(nbins)+hData.GetBinWidth(nbins)
    if ytitle_arg:
        ytitle = ytitle_arg
    else:
        ytitle = 'Events / {:.0f} GeV'.format((xmax-xmin)/(nbins))
    if xtitle_arg:
        xtitle = xtitle_arg
    else:
        xtitle = ''

    if xmin_arg:
        xmin = xmin_arg
    else:
        xmin = hData.GetBinLowEdge(1)
    if xmax_arg:
        xmax = xmax_arg
    else:
        xmax = hData.GetBinLowEdge(nbins)+hData.GetBinWidth(nbins)

    if ymax_arg:
        ymax = ymax_arg
    else:
        all_histos = [h for h in hBkg.values()]
        all_histos.append(hData)
        if dictSig:
            for n, sig in dictSig.items():
                all_histos.append(sig[0])
        all_histos.append(hTot)
        ymax = 1.6 * np.max([h.GetBinContent(i)+h.GetBinError(i)
                             for i in range(1, h.GetNbinsX()+1)
                             for h in all_histos])
    if ymin_arg:
        ymin = ymin_arg
    else:
        ymin = 0
    for b in bkg_name:
        hBkg[b].GetXaxis().SetRangeUser(xmin, xmax)
        hBkg[b].GetXaxis().SetRangeUser(xmax, xmax)
    hData.GetXaxis().SetRangeUser(xmin, xmax)
    hTot.GetXaxis().SetRangeUser(xmin, xmax)

    cwidth, chigh = int(1000*can_scale), int(800*can_scale)
    if plot_ratio:
        cwidth, chigh = int(900*can_scale), int(800*can_scale)
    if can_ratio:
        cwidth, chigh = int(cwidth/can_ratio), chigh
    if canvas:
        canv = canvas
        canv.SetWindowSize(cwidth, chigh)
    else:
        canv = ROOT.TCanvas(plot_name, plot_name, cwidth, chigh)
    ROOT.SetOwnership(canv, False)
    canv.SetTitle('')

    if plot_ratio:
        padhigh = ROOT.TPad('padhigh', 'padhigh', 0., 0.3, 1., 1.0, 0, 0, 0)
        padlow = ROOT.TPad('padlow', 'padlow', 0., 0.0, 1., 0.3, 0, 0, 0)
        padhigh.Draw()
        padhigh.cd()
        padhigh.SetTopMargin(0.08)
        padhigh.SetBottomMargin(0.0)
        padhigh.SetLeftMargin(0.12)
        padhigh.SetRightMargin(0.05)
        padhigh.SetFrameBorderMode(0)
        canv.cd()
        padlow.Draw()
        padlow.cd()
        padlow.SetTopMargin(0.0)
        padlow.SetBottomMargin(0.45)
        padlow.SetLeftMargin(0.12)
        padlow.SetRightMargin(0.05)
        padlow.SetFrameBorderMode(0)
    else:
        padhigh = canv

    padhigh.cd()
    if plot_ratio:
        x1, y1, x2, y2, textsize = 0.61, 0.3, 0.92, 0.90, 0.045
        if leg_ncols==2:
            x1, y1, x2, y2 = 0.43, 0.55, 0.94, 0.90
    else:
        x1, y1, x2, y2, textsize = 0.61, 0.5, 0.98, 0.93, 0.038
        if leg_ncols==2:
            x1, y1, x2, y2, textsize = 0.48, 0.65, 0.94, 0.93, 0.034
    if leg_pos:
        x1, y1, x2, y2 = leg_pos
    if leg_put_nevts:
        textsize = 0.03
    if leg_textsize:
        textsize = leg_textsize

    def make_leg_name(histo, name):
        if leg_put_nevts:
            Etot = ROOT.Double(0.0)
            Ntot = histo.IntegralAndError(-999, 999, Etot)
            if name in ('Data', 'data', 'DATA'):
                return '{} ({:.0f})'.format(name, Ntot)
            else:
                return '{} ({:.0f} #pm {:.0f})'.format(name, Ntot, Etot)
        else:
            return name

    leg = ROOT.TLegend(x1, y1, x2, y2)
    ROOT.SetOwnership(leg, False)
    leg.SetNColumns(leg_ncols)
    leg.SetTextFont(42)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetLineColor(0)
    leg.SetTextSize(textsize)
    leg.AddEntry(hData, make_leg_name(hData, 'Data') , 'lp')
    for b in bkg_name:
        leg.AddEntry(hBkg[b], make_leg_name(hBkg[b], bkg_legname[b]), 'f')
    if dictSig:
        for n, sig in dictSig.items():
            leg.AddEntry(h, make_leg_name(h, legName), 'l')
            h, color, norm, legName = sig
    leg.AddEntry(hTot, make_leg_name(hTot, unc_leg), 'f')


    hData.SetMarkerStyle(20)
    if plot_ratio:
        hData.SetMarkerSize(1.7*can_scale)
    else:
        hData.SetMarkerSize(2.0*can_scale)
    if m_size:
        hData.SetMarkerSize(m_size)
    hData.SetLineWidth(3)
    hTot.SetFillColorAlpha(1, error_alpha)
    hTot.SetFillStyle(error_fill)
    hTot.SetMarkerSize(0)
    hTot.SetMaximum(ymax)
    hTot.SetMinimum(ymin)
    if plot_ratio:
        hTot.GetXaxis().SetLabelSize(0.0)
    else:
        hTot.GetXaxis().SetLabelSize(0.045)
    if plot_ratio:
        hTot.GetYaxis().SetTitleOffset(1.0)
    else:
        hTot.GetYaxis().SetTitleOffset(1.4)
    if plot_ratio:
        hTot.GetYaxis().SetTitleSize(0.055)
    else:
        hTot.GetYaxis().SetTitleSize(0.045)
    hTot.GetYaxis().SetLabelSize(0.045)
    if xticksInt:
        hTot.GetXaxis().SetNdivisions(hTot.GetNbinsX(), 0, 0, True)
    hTot.Draw('E2')
    hstack.Draw('hist same')
    hData.Draw('Esame')
    hTot.Draw('E2same')
    leg.Draw('same')
    if dictSig:
        for n, sig in dictSig.items():
            sig[0].Draw("hist same")
    hData.Draw("esame")
    if bin_label:
        for i, r in enumerate(bin_label):
            hTot.GetXaxis().SetBinLabel(i+1, bin_label[i])
        hTot.GetXaxis().SetLabelSize(0.0)
    hTot.GetXaxis().SetTitle(xtitle)
    hTot.GetYaxis().SetTitle(ytitle)

    if is_logy:
        padhigh.SetLogy()
        if ymax_arg:
            hTot.SetMaximum(ymax)
        else:
            hTot.SetMaximum(ymax*20)
        hTot.SetMinimum(0.02)

    x0, y0, dy, txt_size = 0, 0, 0, 0
    if plot_ratio:
        x0, y0, dy, txt_size = 0.15, 0.84, 0.07, 0.052
    else:
        x0, y0, dy, txt_size = 0.19, 0.87, 0.06, 0.043
    if atlas_label == 'ATLAS':
        ATLASLabel(x0, y0, '', plot_ratio, can_ratio)
    else:
        ATLASLabel(x0, y0, atlas_label, plot_ratio, can_ratio)
    stampText('#sqrt{s} = 13 TeV, '+'{:.1f} '.format(lumi)+'fb^{-1}', x0, y0-dy, txt_size)
    if plot_labels:
        for i, l in enumerate(plot_labels):
            stampText(l, x0, y0-(i+2)*dy, txt_size)

    ROOT.gPad.RedrawAxis()

    if plot_ratio:

        if (ratio_type == 'SoverB' or ratio_type == 'signif'):
            if not dictSig:
                raise NameError('ratio_type \'SoverB\' is not supported when no signal is specified')
            hmc_err = hTot.Clone("hmc_err")
            ROOT.SetOwnership(hmc_err, False)
            hsig, color, norm, legName = list(dictSig.values())[0]
            hmc_err.SetFillStyle(0)
            hmc_err.SetLineColor(color)
            hmc_err.SetLineWidth(3)
            for ii in range(1, hmc_err.GetNbinsX()+1):
                s, b, berr = hsig.GetBinContent(ii), hmc_err.GetBinContent(ii), hmc_err.GetBinError(ii)
                berr2 = berr**2
                if ratio_type == 'SoverB':
                    hmc_err.GetYaxis().SetTitle('S/#sqrt{B}')
                    SoverB = s/np.sqrt(b+berr**2)
                    if np.isnan(SoverB) or np.isinf(SoverB):
                        hmc_err.SetBinContent(ii, 0.0)
                    else:
                        hmc_err.SetBinContent(ii, SoverB)
                elif ratio_type == 'signif':
                    hmc_err.GetYaxis().SetTitle('Z_{A}')
                    try:
                        term1 = (s+b)*np.log((s+b)*(b+berr2)/(b**2+(s+b)*berr2))
                    except ZeroDivisionError:
                        term1 = 0
                    try:
                        term2 = b**2/berr2*np.log(1+berr2*s/b/(b+berr2))
                    except ZeroDivisionError:
                        term2 = 0
                    sig2  = 2*(term1-term2)
                    if np.isnan(sig2) or np.isinf(sig2) or sig2<0:
                        hmc_err.SetBinContent(ii, 0.0)
                    else:
                        hmc_err.SetBinContent(ii, np.sqrt(sig2))
                else:
                    err = 'ratio_type is only \'SoverB\', \'ratio\' or \'signif\', but not \'{}\''.format(ratio_type)
                    raise NameError(err)

                hmc_err.SetMinimum(0.0)
                hmc_err.SetMaximum(1.5)
                cline = ROOT.TF1('cline', '3', -100, 5000)
                ROOT.SetOwnership(cline, False)
                cline.SetLineWidth(1)

        elif ratio_type == 'ratio':
            hdataovermc = hData.Clone()
            ROOT.SetOwnership(hdataovermc, False)
            hdataovermc.Divide(hTot)
            hdataovermc = remove_0entry_data(hdataovermc, 0.01)
            hmc_err = hTot.Clone("hmc_err")
            ROOT.SetOwnership(hmc_err, False)
            for ii in range(1, hmc_err.GetNbinsX()+1):
                if hmc_err.GetBinContent(ii) < 0.001:
                    hmc_err.SetBinError(ii, 0.)
                    hdataovermc.SetBinError(ii, 0)
                else:
                    hmc_err.SetBinContent(ii, 1.0)
                    hmc_err.SetBinError(ii, hmc_err.GetBinError(ii)/hTot.GetBinContent(ii))
                    hdataovermc.SetBinError(ii, hData.GetBinError(ii)/hTot.GetBinContent(ii))
            hmc_err.SetFillStyle(error_fill)
            hTot.SetFillColorAlpha(1, error_alpha)
            hdataovermc.SetMarkerStyle(20)
            hdataovermc.SetLineWidth(2)
            hmc_err.SetMinimum(0.0)
            hmc_err.SetMaximum(2.0)
            hmc_err.GetYaxis().SetTitle("Data / Pred.")
            cline = ROOT.TF1('cline', '1', -100, 5000)
            ROOT.SetOwnership(cline, False)
            cline.SetLineWidth(1)

        else:
            err = 'ratio_type is only \'SoverB\', \'ratio\' or \'signif\', but not \'{}\''.format(ratio_type)
            raise NameError(err)

        padlow.cd()
        if r_ymin:
            hmc_err.SetMinimum(r_ymin)
        if r_ymax:
            hmc_err.SetMaximum(r_ymax)
        hmc_err.GetXaxis().SetTitleSize(0.15)
        hmc_err.GetXaxis().SetTitleOffset(0.9)
        hmc_err.GetXaxis().SetLabelSize(0.12)
        hmc_err.GetYaxis().SetTitleOffset(0.4)
        hmc_err.GetYaxis().SetTitleSize(0.12)
        if xlabel_size:
            hmc_err.GetXaxis().SetLabelSize(xlabel_size)
        if xlabel_offset:
            hmc_err.GetXaxis().SetLabelOffset(xlabel_offset)
        hmc_err.GetYaxis().SetLabelSize(0.12)
        hmc_err.GetYaxis().SetNdivisions(504)
        if ratio_type == 'ratio':
            hmc_err.Draw('E2')
            hdataovermc.Draw('E0 same')
        elif ratio_type == 'SoverB' or ratio_type == 'signif':
            hmc_err.Draw('hist')
        else:
            err = 'ratio_type is only \'SoverB\' or \'ratio\', but not \'{}\''.format(ratio_type)
            raise NameError(err)
        cline.Draw('same')


    if plotdir:
        import os
        full_path_plot = plotdir+'/'+plot_name
        if not os.path.isdir(plotdir):
            os.makedirs(plotdir)
    else:
        full_path_plot = plot_name

    canv.Update()
    canv.SaveAs(full_path_plot+'_{}.pdf' .format(atlas_label))
    canv.SaveAs(full_path_plot+'_{}.png' .format(atlas_label))
    canv.SaveAs(full_path_plot+'_{}.root'.format(atlas_label))
    return canv
