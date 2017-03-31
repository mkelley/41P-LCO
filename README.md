# 41P-LCO
Data pipeline for 41P/Tuttle-Giacobini-Kresak LCO outburst project.

# Project overview
Comet activity varies as different parts of the nucleus are illuminated by sun light, tracing het- erogeneities. Comet 41P will pass Earth very closely this Spring, and is therefore a target of high scientific interest. LCO provides the unique opportunity to acquire high-quality, high-cadence pho- tometry of comets, obtained at uniform-cadence with the same instrument. We request 144 h of observing time on the LCOs 1-m telescopes to reveal rotational, seasonal, and evolutionary changes in the comet’s activity. We will acquire 2-day sampling of its dust (r’) and gas (g’), and will trigger a follow-up program when we detect evidence of outbursts. This dataset will allow us to develop, test, and refine algorithms for the automated detection of outbursts in LSST data, which will be crucial to interpret comet light curves amid the vast amount of data that the survey will produce. 41P’s apparition provides a proxy for a broad range of observing conditions. We will use this program to prototype an early warning system to detect small outbursts, and to activate rapid follow-up observations.

## Team
Survey Team: D. Bodewits, Michael S. P. Kelley, Matthew M. Knight, Tony L. Farnham, Silvia Protopapa, Jian-Yang Li, Tim Lister  
Follow up Team: Quan-Zhi Ye, Henry Hsieh, Zhong-Yi Lin, Emmanuel Jehin, Colin Snodgrass, Alan Fitzsimmons

# Requirements
* Python3
* [astropy](https://www.astropy.org)
* [requests](http://docs.python-requests.org/en/master/)
* [astroquery](https://github.com/astropy/astroquery) >=0.3.5
* [callhorizons](https://github.com/mommermi/callhorizons)
* [matplotlib](https://www.matplotlib.org)

# Install
The recommended installation is through a virtualenv installation.  See the wiki [for details](https://github.com/mkelley/41P-LCO/wiki/Pipeline-installation).

# Pipeline overview
## Inputs
e11 and e91 data

## Outputs
Comet properites
