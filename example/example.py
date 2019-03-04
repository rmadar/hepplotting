import ROOT
import hepplotting as plt

counter = -1

bkg_name = ['bkg1','bkg2','bkg3']

bkg_legname = {
    'bkg1': '#chi#bar{#chi} #rightarrow MM',
    'bkg2': 'pp #rightarrow #psi#bar{#psi}'  ,
    'bkg3': 'H^{+}_{2} with #alpha=1/137'   ,
}

bkg_color = {
    'bkg1' : 868,
    'bkg2' : 867,
    'bkg3' : 866,
}


# ======= CREATE HISTOGRAMS TO RUN THE EXAMPLE ======
def get_random_histo(name, N):
    global counter
    h=ROOT.TH1F(name,name,30,-5,10)
    h.FillRandom('gaus', N)
    h.SetName(h.GetName()+'_tmp{}'.format(counter))
    counter+=1
    return h

n_evts = [246, 251, 252]

dictBkg = {b: [get_random_histo(b, n), bkg_color[b], bkg_legname[b]] for b, n in zip(bkg_name, n_evts)}
dictSig = {'s': [get_random_histo('s', 52), ROOT.kRed+1, 20, 'M_{Madaron}=1 MeV']}
hData   = plt.sum_histograms([get_random_histo('Data', 250) for i in range(0,3)])
hTot    = plt.sum_histograms( [v[0] for v in dictBkg.values()] )
# ==================================================


plt.make_nice_canvas(dictBkg,hTot,hData,plot_name='Example_plot', dictSig=dictSig,
                     ytitle='Probability Density Function',
                     xtitle='Random variable', plot_ratio=True,
                     ymax=300, ratio_type='signif',
                     leg_ncols=1, leg_put_nevts=True, leg_textsize=0.036)
