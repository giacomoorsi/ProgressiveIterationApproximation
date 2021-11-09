
# Libreria di funzioni di approssimazione e interpolazione PIA su B-Spline in MATLAB

Questa cartella contiene il codice che ho implementato per la mia tesi di laurea in Informatica per il Management all'Università di Bologna, conseguita nel 2021. 

**Titolo tesi**: Progressive Iteration Approximation: una classe di metodi per interpolazione e approssimazione spline  
**Relatore**: prof. Giulio Casciola  
**Presentata da**: Giacomo Orsi

La tesi è accessibile nel file `tesi.pdf`

## Abstract
Lo scopo di questa tesi è introdurre alcuni metodi iterativo-geometrici di interpolazione e approssimazione per curve spline e analizzarne le differenze. In particolare, si analizza una classe di metodi denominati "Progressive Iteration Approximation". Questi metodi permettono di determinare i punti di controllo delle curve interpolanti o approssimanti in modo iterativo con il vantaggio di gestire insiemi di dati molto grandi in modo efficiente. Inoltre, si presenta una sintesi di alcuni risultati numerici ottenuti in ambiente MATLAB a seguito dello studio, dell'implementazione e della sperimentazione degli algoritmi analizzati.


## Codice 
Questa cartella contiene il codice relativo all'implementazione per curve B-Spline in MATLAB degli algoritmi
- PIA (Progressive Iteration Approximation) (https://doi.org/10.1016/j.cad.2013.08.012)
- WPIA (Weighted Progressive Iteration Approximation) (https://doi.org/10.1016/j.cagd.2009.11.001)
- PIA Adaptive (https://doi.org/10.1016/j.cagd.2012.03.005)
- LSPIA (Least Squares Progressive Iteration Approximation) (https://doi.org/10.1016/j.cad.2013.08.012)
- LSPIA Progressive (https://doi.org/10.1016/j.cad.2013.08.012)


PIA e WPIA sono entrambi contenuti all'interno della cartella `pia_bspline`. 

Ogni cartella contiene il file `<nomeAlgoritmo>_body.m` e  `<nomeAlgoritmo>_example.m`

All'inizio del file `body` sono mostrati tutti i parametri da impostare nell'ambiente  prima di eseguirlo. 

Nei file `example` sono mostrati degli esempi di parametri per eseguire il `body`. 

Nella cartella `svg` sono presenti diverse immagini svg composte da solo una *path*, e quindi sono approssimabili/interpolabili con una sola curva. 

Nella cartella `librerie` sono presenti le librerie utilizzate. 

Per utilizzare gli script, si suggerisce di aprire tutta la cartella in MATLAB e eseguire lo script `load_libraries.m` che aggiunge alla working path tutte le librerie e le immagini svg. 

Nella cartella `figure_tesi` è presente un file `.txt` con le istruzioni per replicare le figure mostrate nel documento di tesi. 