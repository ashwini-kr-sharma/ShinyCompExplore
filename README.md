# ShinyCompExplore

This is a Shiny app to run supervised and unsupervised cell type deconvolution methods. It aids in a user-friendly analysis and visualization of various deconvolution algorithms.

### Connecting `ShinyCompExplore` with the Docker image `mrna_meth_decon` and starting the app

[`mrna_meth_decon`](https://hub.docker.com/r/ashwinikrsharma/mrna_meth_decon) is a Docker image which encompasses the following deconvolution algorithms and Shiny along with their dependencies.

immunedeconv\
estimate\
CellMix\
DeconRNASeq\

EpiDISH\
EDec\
medepir\

ica\
deconica\
fastICA\
NMF\

```
cd ~/

docker pull ashwinikrsharma/mrna_meth_decon

git clone https://github.com/ashwini-kr-sharma/ShinyCompExplore.git

docker run --rm -it -p 3838:3838 -v ~/ShinyCompExplore:/srv/shiny-server/ ashwinikrsharma/mrna_meth_decon

# Wait for few seconds for the shiny app to start

```
