# High quality HEP plotting

This code allows to produce, in a very flexible way, publication quality HEP distributions
based on ROOT histograms.

## 1. Installation

```
git clone https://github.com/rmadar/hepplotting
cd hepplotting
pip install -e . --user
```

## 2. Usage

### 2.1 What do you need

1. Histograms for every processes and data, and total histogram (for the uncertainty)
2. Legend names and color definitions for every background
3. Possibly every systematic variations to be passed through the total histogram
(after appropriate combination across systematics)

### 2.2 Simplest call

The simplest use you can do is to call the following function, once you get your
histograms (with their color):
```
import hepplotting as plt
plt.make_nice_canvas(dictBkg,hTot,hData,plot_name='myplot')
```
where:
  + `dictBkg` is a dictionnary made of background name and the corresponding histogram `{bname:[TH1F,color,leg_name]}`
  + `hTot` is the total histogram with possibly larger uncertainty (to account for systematics)
  + `hData` is the data histograms

A concrete example can be found in [SimpleExample.py](SimpleExample.py) which produces this plot:
![distribution](plots/Example_plot_Internal_reduced.png)



## Technical comments

### Dependencies

  + ROOT
  

### Known issues

The notebook example doesn't seem to work well and all setup. Several minor features are not working properly
  + the magic command `autoreload` doesn't work (nothing to do with this pacakge)
  + the displayed canvas is not correct (should be connected to ROOT-related tool versions). This affect the 
  use of the tool within a notebooks, only.
