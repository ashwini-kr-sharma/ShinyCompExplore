# ShinyCompExplore

This is a Shiny app to run supervised and unsupervised cell type deconvolution methods. It aids in a user-friendly analysis and visualization of various deconvolution algorithms.

### Connecting `ShinyCompExplore` with the Docker image `mrna_meth_decon` and starting the app

[`mrna_meth_decon`](https://hub.docker.com/r/ashwinikrsharma/mrna_meth_decon) is a Docker image which encompasses the following deconvolution algorithms and Shiny along with their dependencies.

(1) immunedeconv\
(2) estimate\
(3) CellMix\
(4) DeconRNASeq

(5) EpiDISH\
(6) EDec\
(7) medepir

(8) ica\
(9) deconica\
(10) fastICA\
(11) NMF

```
cd ~/

docker pull ashwinikrsharma/mrna_meth_decon:latest

git clone https://github.com/ashwini-kr-sharma/ShinyCompExplore.git

docker run --rm -it -p 3838:3838 -v ~/ShinyCompExplore:/srv/shiny-server/ ashwinikrsharma/mrna_meth_decon

# Or, for deployment in a Virtual Machine with a public IP
# docker run -d -p 3838:3838 -v ~/ShinyCompExplore:/srv/shiny-server/ ashwinikrsharma/mrna_meth_decon

# Wait for few seconds for the shiny app to start

```

Once you see a message `Starting listener on http://[::]:3838`, go to your browser (Chrome, Firefox etc) and type in the address bar `localhost:3838`, wait for some time for the page to load. Now, you should now be able to see the `ShinyCompExplore` app.

__Demo data__ is made available inside `~/ShinyCompExplore/data/demo`, that can be used to test the app.

Raise a [Github issue](https://github.com/ashwini-kr-sharma/ShinyCompExplore/issues) in case you face any problems !!
