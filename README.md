# EDiSC
This repository contains the R scripts and data files used to produce the results reported in Zafar and Nicholls (2023) "An Embedded Diachronic Sense Change Model with a Case Study from Ancient Greek".

**Data extraction and snippets:** The ancient Greek data used to extract snippets for the target words "kosmos", "mus" and "harmonia" is freely available online. Links included in the file _"Greek data extraction.R"_ point to the data sources, and the code may be used to extract the snippets. Alternatively, the files _"kosmos_snippets.RData"_, _"mus_snippets.RData"_ and _"harmonia_snippets.RData"_ contain the extracted snippets ready to use. 

The English data used to extract snippets for the target word "bank" needs to be purchased from https://www.english-corpora.org/coha/. (Note that the COHA data used in this paper was the version as at 3 June 2020.) The file _"COHA data extraction.R"_ may be used to extract the snippets once the data has been obtained. With permission from the data publisher, the extracted snippets are included in the file _"bank_snippets.words.RData"_, as well as our manual sense annotation for these snippets, ready to use.

**Embeddings:** The files _"Embeddings - Greek data.R"_ and _"Embeddings - COHA.R"_ can be used to generate the GloVe embedding vectors for the ancient Greek and COHA data respectively. Alternatively, _"word vectors"_ contains the embedding vectors for all context words (after data filtering) for the four target words ready to use.

**Models and samplers:** R scripts used to fit the models using HMC and MALA (not Stan), and using NUTS and ADVI (Stan). 

**Scalability analysis:** R scripts used to simulate the data and perform the runs used in Figure 9 in the paper.

**figures and tables.R:** Contains the code used to produce all the figures, tables and results reported in the paper.

