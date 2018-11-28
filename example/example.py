import ROOT
import hepplotting as plt

counter = -1

bkg_name = ['bkg1','bkg2','bkg3']

bkg_legname = {
    'bkg1': 'Bkg_{1} (m_{#chi}=10 GeV)',
    'bkg2': 'Bkg^{2}'                  ,
    'bkg3': 'Bkg_{3} #alpha=1/137'     ,
}

bkg_color = {
    'bkg1' : 868,
    'bkg2' : 867,
    'bkg3' : 866,
}


# ======= CREATE HISTOGRAMS TO RUN THE EXAMPLE ======
def get_random_histo(name):
    global counter
    h=ROOT.TH1F(name,name,30,-5,10)
    h.FillRandom('gaus',250)
    h.SetName(h.GetName()+'_tmp{}'.format(counter))
    counter+=1
    return h

dictBkg = {b:[get_random_histo(b),bkg_color[b],bkg_legname[b]] for b in bkg_name}
dictSig = {'s': [get_random_histo('s'), ROOT.kRed+1, 20, 'Madaron (m_{#xi}=1 MeV)']}
hData   = plt.sum_histograms([get_random_histo('Data') for i in range(0,3)])
hTot    = plt.sum_histograms( [v[0] for v in dictBkg.values()] )
# ==================================================


plt.make_nice_canvas(dictBkg,hTot,hData,plot_name='Example_plot', dictSig=dictSig,
                     ytitle='Probability Density Function',
                     xtitle='Random variable',plot_ratio=True,
                     ymax=300, ratio_type='SoverB')
