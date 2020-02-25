---
Enrico Verdolotti, 
January/February 2020, 
Final project for parallel and distributed computing course

# Molecular Dynamics with CUDA
## meta-project

This project will be based on a previous work of Emily Crabb, publicly
available at: <https://github.com/ejc44/MD>. The main objective is to successfully
implement the same algorithms and specifications, harnessing GPU
multi-threading capability and in particular, exploiting CUDA toolkit.
Moreover will be done some review of the original code and introduced
new stuffs as described below.

The envisioned steps are:

1.  Review and rewrite **serial** version;

2.  Write **GPU multi-thread** version with **CUDA**;

3.  Review and rewrite **CPU multi-thread** version;

4.  **Benchmark** and comparisons;

5.  New **animated visualization**.

Two years are passed since the original work so the revision of
**serial** version will be done mainly for better understanding of
algorithm, data structures involved, entire workflow and, possibly for
old-code revision using new constructs that maybe was not available in
the previous version of Julia, taking the opportunity for refactoring
and finding bottlenecks (e.g. the \"find forces\" function, re-allocate
for each time step the entire forces \"array\" and maybe this can be
avoided.)

In the original work, the **GPU version** implementation was left
incomplete due to the complexity of main functions which contains
several conditionals check. That make hard to write GPU code using just
built-in functions and CuArrays. The idea for addressing this problem is
to write an ad hoc kernel function capable to run on GPU along with the
classic CuArrays usage.

It would be nice to have a **CPU multi-thread** version more similar to
GPU version for comparisons and even if the common workflow for
developing multi-threaded code is to write first CPU multi-thread
version and then the GPU one, in this case will be done in reverse due
to two reasons: First, the multi thread CPU code is already available
and working and, second, to adapt a kernel function to CPU or GPU from
the respective counterpart is relatively easy.

As ever, will be mandatory to test if a real speedup has been achieved.
For this reason,\
**BenchmarkTools** package will be used to compare both current results
and original ones.

Lastly, a more expressive **visualization** of particle motion will be
explored using classic Plot library or Makie for 3D animations.

In conclusion some of the original project, such as distributed
computing, part will not be considered even if these was written with
`@parallel` macro that is no longer supported and this would need a code
rewriting. That choice has been done due to the low score obtained for
distributed version in original code and lack of time but that doesn't
mean it cannot be resumed in future.

> *"We adore chaos because we love to produce order."* -- *M.C. Escher*

---
--- 
---
Enrico Verdolotti, 
Gennaio/Febbraio 2020,
Progetto finale per il corso di calcolo parallelo e distribuito 

# Molecular Dynamics with CUDA
## meta-progetto

Questo progetto si baserà su un precdente lavoro di Emily Crabb,
diponibile all'indirizzo: <https://github.com/ejc44/MD>. 
L'obiettivo principale è implementare con successo gli stessi algoritmi e specifiche
sfruttando le capacità di multi-threading della GPU ed in particolare con il toolkit CUDA.
Inoltre verrà effettuata una revisione del codice originale ed introdotte delle novità
come descritto a seguire.

Le tappe previste sono:

1.  Rivedere e riscrivere la versione **seriale**;

2.  Scrivere la versione **GPU multi-thread** con **CUDA**;

3.  Rivedere e riscrivere la versione **CPU multi-thread**;

4.  **Benchmark** e comparazione delle prestazioni;

5.  Nuove **visualizzazioni animate**.


Sono trascorsi due anni dal lavoro originale, perciò la revisione
della versione **seriale** sarà effettuata principalmente per una
migliore comprensione degli algoritmi, strutture dati utilizzate,
intero flusso di lavoro ed eventualmente revisione del vecchio codice,
utilizzando nuovi costrutti che forse non erano disponibili nelle precedenti
versioni di Julia, cogliendo l'opportunità per eseguire il refactoring
e scoprire colli di bottiglia (e.g. la funzione \"find_forces\" ri-alloca
ad ogni istante della simulazione l'intero \"array\" delle forze e forse
questo può essere evitato.)


Nel lavoro originale, l'implementazione della **versione GPU** è stata
lasciata incompleta a causa della complessita della funzione principale che
contien diversi controlli condizionali. Ciò rende difficile la scrittura di
codice GPU utilizzando solo le funzioni integrate ed i CuArrays. L'idea per
risolvere questo problema è di scrivere una funzione kernel apposita in grado 
di girare sulla GPU assieme al classico uso dei CuArrays.


Non sarebbe male avere una versione **CPU multi-thread** più simile alla versione
GPU per confronto, ed anche se normalmente per sviluppare codice multi-threaded si 
scrive prima la versione per CPU e poi quella per GPU in questo caso verrà fatto 
il contrario per due ragioni: Primo, il codice multi thread è già disponibile e funzionante
e, Secondo, adattare una funzione kernel alla GPU od alla CPU dalla rispettiva controparte
è relativamente facile.


Come sempre sarà obbligatorio testare se si è ottenuto un miglioramento effettivo.
A questo scopo, verrà utilizzato il package **BenchmarkTools** per confrontare
sia i risultati correnti che quelli originali.

Infine, verrà esplorata una nuova possibilità di **visualizzazione**
del moto delle particelle usando la classica libreria Plot oppure Makie 
per le animazioni 3D.

Concludendo, alcune parti del progetto originale, come la versione
relativa al calcolo distribuito, non saranno cosiderate anche se sono
state scritte con la macro `@parallel` che non è più supportata e perciò
avrebbe bisogno di una riscrittura. Questa scelta viene fatta a causa del
risultato negativo dei benchmark relativi alla versione distribuita del precedente 
lavoro e per mancanza di tempo ma ciò non vuol dire che non possa essere 
ripresa in futuro.

> *"We adore chaos because we love to produce order."* -- *M.C. Escher*

---
---

# Report

### 1. Rivedere e riscrivere la versione seriale
In scrittura...

### 2. 

### 3.

### 4.

### 5.

