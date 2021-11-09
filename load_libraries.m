% Determino la cartella nella quale è situato il file
folder = fileparts(which('load_libraries.m')); 
% Aggiungo quella cartella e tutte le sottocartelle
addpath(genpath(folder));