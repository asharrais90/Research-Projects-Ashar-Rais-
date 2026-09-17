**Contents:**

**1.** **Introduction**

	1.1.     Pluripotency stages and chromatin dynamics
	
	1.2.    Tip-Enhanced Raman Spectroscopy (TERS) and Raman Dots (R-dots)

**2.** **Major objectives**

**3.** **Methodology**

	**3.1.**     **Cell culture for pluripotent stem cells**
	
	**3.2.**     **Nuclear isolation**
	
	**3.3.**     **Fixation & Tokuyasu cryo-sectioning**
	
	**3.4.**     **Preparation of R-dots & TERS**

**4.** **Discussion**

	**4.1.**     **Comparison with other imaging systems**
	
	**4.2.**     **Limitations of ct-TERS**
	
	**4.3.**     **Future Prospects**

  

**1.** **Background & Motivation**

**1.1.**      **Pluripotency stages & chromatin dynamics**

Studying the development of our body has been one of the oldest endeavours in human history. Different models have been proposed over the centuries to try and explain this mysterious process. With the advent of modern scientific techniques, the details of this process are being discovered at an unprecedented rate and we now know enough to even create systems in a lab that faithfully recapitulate some aspects of this process**1**. However, what we don’t understand still dominates over what we have come to understand and to turn this around, great amount of work is necessary. The project proposed here is part of this combined effort to understand how we are created.

Humans, like almost all other organisms, start their developmental journey as a single cell called as the Zygote. This cell has the potential to form every cell required for the developmental process (embryo plus extra-embryonic tissue like the placenta) and is considered to be ‘totipotent’ (gr. Toti- total & potent- potential). As this single cell divides to create the required number of cells, the potential of each cell decreases and is traded off for the ability to perform a specific cellular function. So, after totipotency comes pluripotency and cells with this level of potential now have the capacity to form embryonic tissue only**2**. When looking at the developmental stages, the cells having this potential are present from the blastocyst stage (specifically the cells of the inner cell mass or ICM) to just before the advent of somitogenesis (around 2 weeks in human development)**3**.  This stage is of particular importance because it encompasses events like implantation into the uterus and gastrulation (formation of the 3 germ layers), where many things can go wrong leading to failed pregnancies and congenital defects**4**. Therefore, it is of great significance that we understand the state of pluripotency and how the development of embryo at this stage proceeds without any hiccups. Another significance of the pluripotent state is its role in regenerative medicine, with certain protocols existing to reprogram differentiated cells (like fibroblasts) back to pluripotency. This has various benefits including developmental monitoring, disease modelling, replacement of lost/faulty tissues, etc**5,6**. Considering all of these points, the objective of understanding pluripotency becomes like a crucial duty for the scientific community and this project will aid in furthering us towards achieving this objective.

Over the years, the data collected from a cohort of experiments points to the fact that pluripotency itself is divisible into sub-states, each having their own subtle characteristics and importance for specific events. In 2009, J. Nichols & A. Smith proposed two sub-states: naïve and primed pluripotency**7**. Former is the state of the pluripotent cells before implantation, while latter is post-implantation status of these cells. Morphologically, the colonies formed by each sub-type are different, with naïve colonies having a smooth dome-shaped appearance and primed colonies, a flatter structure with rougher boundaries.  This sub-division is of particular relevance for culturing of pluripotent stem cells in-vitro and making cellular aggregates for modelling embryonic development. For some reason, when cells are taken from pre-implantation human embryos and cultured in-vitro, they demonstrate the characteristics of primed stem cells instead of naïve stem cells. This is important because the culturing conditions is different for both of these cell types and this difference lies mainly in the types of signalling modulators used in the medium**8**. Following years of study of these two sub-types, it became apparent that another sub-type (that is distinct enough from the continuum-like nature of other 2 sub-types) is present between them- the formative pluripotency. Hypothesized again by A. Smith in 2017**9** the cellular and sub-cellular features of the formative state have since been elucidated by many follow up studies**10–12**, and this state represents the transition between the naïve and primed states, making it crucial for creation of certain embryonic models.

At the sub-cellular scale, many differences are observed between these three sub-types**9,13**. Some of these include:

·       DNA methylation status: As cells progress from naïve to formative, then primed state, an increase in global DNA methylation of the genome has been observed**9**;

·       X-inactivation: Naïve (female) cells have 2 active X-chromosomes while primed cells have one inactive and one active. This X-inactivation process apparently happens during the formative state**9**;

·       Differential expression of naïve (Klf2/4, Tbx3, etc.), early post-implantation (Sox3, Otx2, etc.)  and lineage factors (cer1, foxa2, brachyury, etc.)**9**;

·       Specific responses to different signalling axes (like TGF-β and FGF-ERK signalling pathways)**14**;   

·       Chromatin dynamics: (explained below);

·       Germ cell competence: Only the formative state is considered to have the capacity to give rise to primordial germ cells (PGCs), with the naïve state being too early and the prime state being too late**9**;

Chromatin dynamics of pluripotent stem cells is dependent on the type of environmental cues which direct them towards a specific cell lineage/fate (epigenomic modulations)**15**. As already mentioned above, the global genomic DNA methylation changes between the different pluripotency sub-states and this is important to differential gene expression patterns (also mentioned above). In addition to DNA methylation, other types of chromatin modifications play crucial roles in modulating the expression profiles, with histone modifications displaying the most versatile and crucial roles. These modifications, along with the functioning of architectural remodelling proteins (like CTCF, Cohesin, etc.), results in very rapid yet sufficiently distinct organisational patterns of the chromatin (namely TADs and chromatin loops)**16,17**. Of particular interest at this stage of development is a type of chromatin modification observed on the promoters of genes involved in cellular identity functions. Named as the “poised or bivalent domains”, this modification involves the presence of H3K4me3 (an expression activating mark) and H3K27me3 (an expression suppressing mark) and a stalled RNA polymerase II**17**. These poised domains are important as they make it possible for the cell to very rapidly respond to specific developmental cues by activation and suppression of necessary genes (both of which only require the removal of one histone modification). Studying the features of chromatin dynamics and other epigenomic properties holds great promise in advancing our understanding of pluripotent stem cells and will result in better protocols for culturing these sensitive cells in vitro.     

**1.2.**     **Tip-Enhanced Raman Scattering (TERS) & Raman Dots (Rdots)**

Raman spectroscopy is a vibrational spectroscopic technique which probes the sample using light (mostly in the visible and near-IR region) and is based upon the principle of inelastic or raman scattering. To briefly describe the phenomenon of raman scattering: When an electron (in a chemical bond) absorbs a photon, it can either release a photon of the same energy (elastic scattering) or release a photon of different energy (inelastic or raman scattering). For the latter type of scattering, the energy difference between incident and released photons is dependent upon the configuration of vibrational energy states of that molecule. This configuration is unique as it is based on the constituent atoms, properties of the chemical bonds, particular movement of the bond (like stretching, bending, etc.), environment, etc. So based on the profile of raman scattered photons, essentially a “molecular fingerprint” is obtained (called the raman spectrum) and this allows very accurate characterization of a plethora of molecules present in the sample of interest. The major problem with this spectroscopic probing technique is that raman scattering is inherently a very weak phenomenon when compared to elastic scattering as only 1 out 106-108 photons absorbed are scattered in this way**18**. This not only results in issues with sensitivity, but also the time required for the acquisition of sufficient number of photons to produce a reliable spectrum/image. Another issue lies in the spatial resolution offered by the classical raman imaging techniques. Like optical microscopy, the spatial resolution obtained depends upon the wavelength of the probing light and the properties of the imaging system (Abbe’s diffraction limit: spatial resolution ~= λ/2) and is ~300 nm**19** when considering the various parameters. These issues make classical raman spectroscopy not very suitable for the purposes of visualising sub-cellular processes.

To overcome the issues mentioned above, many different methods have been developed to modify different aspects of raman spectroscopy. One of these methods is the Tip Enhanced Raman Spectroscopy (TERS). This technique brings to classical raman system the advantages of Local Surface Plasmon Resonance (LSPR), which is basically the enhancement of electromagnetic fields (or light) around a surface (with specific properties) due the oscillations of groups of electrons present in the surface (and these groups are called plasmons). Due to LSPR, field enhancements of upto 1012 can be observed in very small regions (few nanometres.) around the illuminated surface**20**. In TERS, the surface that is illuminated to get field enhancements is a small tip (nanometre-scale dimensions), attached to a cantilever. This apparatus is scanned over the surface of the sample at a very minute separation, so the amplified raman signal is obtained from very specific spots where the tip is present, greatly reducing the background noise**21**. These factors drastically enhance the quality and quantity of raman signal obtained, allowing for single molecule detection, data acquisition speed in milliseconds and spatial resolution down to 1.7 nm**22**. Additionally, these parameters are tunable by the modulation of many different technical aspects of the system such as the dimensions of the tip, properties of the light used for illumination, the substrate upon which the sample is placed, the detector used for data collection, etc. All of this make TERS a technique with incredible potential for sub-cellular probing studies, and some very interesting uses (like DNA sequencing**23** and imaging of single colon cancer cell surface**24**) have already been demonstrated.

Although one of the main attractions of raman spectroscopy in biological setting is the “probe-less imaging”, it still allows the use of molecules with specific raman spectra as probes. These molecules include alkynes, nitriles or just stable isotope containing species and are named as Raman tags. They have the capacity to greatly increase the utility of raman imaging and numerous applications of raman tags for biological purposes have been developed such as tissue imaging and single cell metabolomic imaging**19**. One type of raman tags, the raman dots (Rdots), created by enclosing many molecules of a raman tag in a polymer-based nanoparticle, have recently garnered a lot of attention**25**. This concentration of many probes in a small volume (~20 nm) drastically increases the brightness of these probes, amplifying the imaging capabilities. A useful property of Rdots is that their excitation and emission spectra are significantly narrower (~50 times) than the probes used for fluorescent microscopy (this technique will be explained in a later section). This allows for multiplexing (which is the simultaneous use for multiple probes for imaging a single sample) of upto 20 channels at once in vitro**25** and 14 channels in live-cells**26** and that is to a much greater degree than fluorescent microscopy (only 4-5 channels**25**).  Bio-functionalization of Rdots is also possible with the use of certain coating materials like PEG (polyethylene glycol), immunoglobulins, amine-modified DNA molecules, etc. and this permits their use for labelling of specific sub-cellular targets**26**.                                      

**2.** **Major objectives:**

·       Development of a novel imaging system, named as cryo-tomographic TERS (ct-TERS), which combines the principles of nuclear cryo-sectioning and TERS.

·       Use of ct-TERS to characterise the chromatin dynamics of the pluripotency sub-states (naïve, formative and primed) which will be presented as a 3-D. This model would include the physical organisation of chromatin fibres and the localisation of transcription factors, epigenetic modulators (like PRC complexes), RNA molecules (long-coding RNA like Xist) and RNA polymerase II at spatial resolutions expected to be around 10-15 nm.

**3.** **Methodology**

The work flow will be divided in two general phases with multiple steps:

·       Isolating nuclei from cells-of-interest, fixing them and slicing them using the Tokuyasu cryo-sectioning method;

·       Introducing Rdot-containing solution to the cryo-sections before using TERS to scan them, obtain the image of each sample. The images from different slides will be combined using specific algorithms to create a 3D model of the chromatin.

**3.1.**      **Cell culture for pluripotent stem cells**

For culturing pluripotent stem cells, the protocol mentioned here will be used**27**. This protocol uses H9 hESC cell lines and employs either chemical cocktail (HDAC + WNT inhibitors) or a 5i/L/A + NK2 transgene method**28** to convert the primed cells back to a naïve state. This conversion takes 2-3 weeks and 10 cells will be collected every 48 hours for ct-TERS, starting just before the first treatment with the conversion cocktail, and the final collection will be when the cells are confirmed to have achieved the naïve state (so an estimated 70-105 samples).

**3.2.**      **Nuclear isolation**

With the recent advent of micro-fluidics based single cell and single-nuclei focussed techniques, more and more of the nuclear isolation methods are being developed and optimised. For this project, individual nuclei will be isolated based on the protocol described by 10X Genomics Single Cell protocol**29**.

**3.3.**      **Fixation & Tokuyasu cryo-sectioning**

The isolated nuclei will be fixed and sliced at very low temperatures (Tokuyasu method) according to the protocol mentioned in this previous study**30**. The fixation step will involve the use of para-formaldehyde (4% for 10 min and 8% for 2h) in HEPES buffer. After fixation, the nuclei will be embedded in sucrose-PBS solution (a cryo-protectant) before plunging them into liquid nitrogen at a very rapid pace. This rapid freezing will prevent the formation of disruptive structures (like ice crystals) and keep the ultrastructure intact. These frozen samples will then be sliced using an ultracryomicrotome, producing nuclear sections ideally ~100 nm in thickness. Once these cryo-sections are produced, further steps can be done at room temperature or even at 37 oC (as done by this group during incubation of cellular cryosection with FISH probes for 40h**30**).  

**3.4.**      **Preparation of R-dots & TERS**

The synthesis of Rdots will be done as mentioned in this study25, and after necessary bio-functionalization steps, the following 10 nuclear components will be targeted:

·       Proteins (use of immunoglobulins): RNA polymerase II, transcription factors (Oct4, Sox2, Klf4, c-Myc & Nanog), epigenetic modulators (PRC-2 complex), Histone marks (H3K4me3 & H3K27me3) and chromatin remodelling factors (CTCF & Cohesin);

·       RNA molecules (use of oligos): long coding RNA (Xist).

Each of the cryo-sections will be incubated overnight, at room temperature with a solution containing all of these Rdots to allow stable labelling.

For TERS, the setup will be done following the protocols of this study for cell surface mapping**24**. This setup includes a raman spectrometer (in backscattering mode) and an AFM (atomic force microscope) apparatus (in non-contact mode). The time for data acquisition for a single cryo-section will range from msecs to secs and the hyperspectral data analysis and unmixing will be done using MATLAB software (Release 2007a). The data from each cryo-section of a particular sample will be combined using z-stacking algorithms used commonly in confocal microscopic analyses, to create a 3-D model of the chromatin.

**4.** **Discussion**

**4.1.**      **Comparison with other imaging systems**

**4.1.1.**      **With Fluorescence microscopy**

Fluorescence microscopy is a type of optical microscopic technique, i.e. it uses light in the visible region (400-700 nm) for visualization purposes. The main working principle of this technique is the use of fluorescent probes to tag desired biological molecules. These fluorescent probes absorb photons of certain wavelengths and after some time, release photons of a higher wavelengths. This property is exploited to get information (concentration, localization, interactions, movement, etc.) only of specific molecules in the cellular crowd. The simplicity and versatility have made fluorescence microscopy among the most popular imaging techniques in biological settings. As fluorescence microscopy uses photons of visible range, the spatial resolution limit for classical systems is ~ 200nm (Abbe’s diffraction limit)**31**. This is the main obstacle in using fluorescent microscopy for the study of sub-cellular processes/structures which are far smaller than this limit of 200 nm. Numerous methods have been developed to successfully surpass this obstacle and techniques which use these methods have been described by the term “super resolution fluorescent microscopy (SRFM)”**32**. Great strides have been made using SRFM, and with one group even achieving angstrom-level of resolution**33**, allowing them to image bases in a DNA strand and study the organisation of CD20 receptor at a scale never achieved before.

When compared to ct-TERS for the task of chromatin architecture modelling, fluorescent microscopy has two significant disadvantages:

·       Multiplexing: As mentioned in a previous section, fluorescent microscopy allows (currently) for only 4-5 simultaneous channel imaging**25**. This, when compared to 14 or more-simultaneous channel imaging for TERS, results in lower capability for studying crucial nuclear components like proteins and RNA molecules.

·       Washout of signal of small molecules (TFs, RNAs, etc.) by the signal of labelled chromatin: For getting a fluorescent signal from chromatin, either probes for DNA (like nucleotide analogues, DAPI, etc.) or probes for histones would have to be used. Both of these methods can result in obtainment of signals from other target molecules becoming problematic as the optical profiles of these kinds of probes is not as narrow as raman tags, and also relative signal strength of commonly used probes are close to each other (unlike the difference between the signal strength of chromatin and Rdots).    

**4.1.2.**      **With Electron microscopy**

In electron microscopy (EM), the wave nature of electrons is for used for imaging purposes**34**. Since electrons have far shorter wavelengths than photons, the resolution power offered by electron microscopy is greater than optical microscopic techniques (ranging from few nm to Å-level) and are regularly employed to image sub-cellular architecture to single protein structure. These kinds of imaging procedures have a particular limitation for biological uses – the imaging must be done in a vacuum (since the electron beam can be disturbed by air particles present between the source and the sample). This effectively puts live-sample imaging out of reach and some even consider live-sample imaging to not be possible. Another limitation of EM is the lack of multiplexing capabilities. While there have been some advances in this area**35**, they are still not close to the multiplexing capabilities offered by fluorescence and raman imaging techniques. In one technique, called as correlative light & electron microscopy (CLEM) fluorescent microscopy is combined with EM to try and exploit the benefits of both the techniques, however many issues arise during this process**36**. So (as of now) for the proposed goal of creating an imaging system for characterizing chromatin architectural changes in pluripotent stages, EM is not suitable mainly because of the lack of multiplexing capabilities required for the study of numerous protein & RNA molecules.

**4.2.**      **Limitations of this system**

Some of the foreseeable limitations of this system include:

·       As it is the case when creating novel systems of any kinds, there are numerous technical challenges which have to be very precisely handled in order to achieve the desired results. These challenges are mainly associated with the interdisciplinary nature of this system, which necessitates the use of methods and protocols which are not usually used together and thus can have a lot of “technical friction” between them. To get through these limitations, a significant number of optimizations and modifications to existing protocols maybe required that can be costly and make the technique less feasible for even the possible uses;

·       Live-cell and time-lapse imaging is not possible with this system as the steps of fixation and cryo-sectioning are integrals part of the operation. This limits the possible uses of this system and also creates a dependency of the results upon the chosen time points for sample collection (so, if a process happens in between the time points, it will escape detection);

·       Even though the protocols for fixation, cryo-sectioning and TERS have been regularly used for different purposes, they are still not completely perfect and can introduce artifacts in the results including damage to the nuclei during isolation, disruption of the chromatin architecture during the cryo-sectioning steps, etc. To address these issues, data from multiple samples will to be compiled to get the safest result possible.

·       The properties of chromatin in pluripotent stem cells are very sensitive to the environmental cues and this can result in variations dependent upon the specific parameters of the stem cell culture.

**4.3.**      **Future prospects**

Understanding epigenomics and especially chromatin dynamics is part of a common trend of this century’s molecular biological endeavours as researchers are now thinking beyond the gene-centric model of molecular biology that dominated the landscape of previous century’s endeavours**37**. As our knowledge of these epigenomic processes expands, we become more comfortable approaching problems that perplexed even the greatest of minds of the past generations (one such problem is how we go from a single cell to a well-formed organism) and a major factor for this is that we understand the need for the development of novel systems for studying these complicated problems. For ct-TERS, regardless of the technical complexities awaiting the development, incredible potential is undoubtedly present. The main strength of this system is the combination of the high spatial resolution and unrivalled multiplexing capabilities, which permits precise modelling of a very complicated environment, and specific probing of many molecules central to controlling the epigenomic processes. This will prove very useful in unravelling the various aspects of pluripotent stem cell biology and provide some form of guidance to many areas of translational medicine, especially regenerative medicine. Additionally, the probing parameters of ct-TERS can be modulated by controlling specific factors of the TERS system and by using different target combinations for the raman probes. Because of these factors, the ct-TERS system is applicable to other important biological questions such as how certain drugs affect the chromatin landscape, what is the role of chromatin dynamics in tumour initiation and progression, etc. All this makes it exciting to see and experience the possible paths opened by ct-TERS, making it worthwhile to invest resources in the development of this imaging system, and in the near future, its improvement.      

**4.4.**      **References**

1.         Zernicka-Goetz, M. The evolution of embryo models. _Nat Methods_ **20**, 1844–1848 (2023).

2.         Du, P. & Wu, J. Hallmarks of totipotent and pluripotent stem cell states. _Cell Stem Cell_ **31**, 312–333 (2024).

3.         Osorno, R. _et al._ The developmental dismantling of pluripotency is reversed by ectopic Oct4 expression. _Development_ **139**, 2288–2298 (2012).

4.         Ferrer-Vaquer, A. & Hadjantonakis, A.-K. Birth Defects Associated with Perturbations in Pre-implantation, Gastrulation & Axis Extension: from Conjoined Twinning to Caudal Dysgenesis. _Wiley Interdiscip Rev Dev Biol_ **2**, 427–442 (2013).

5.         Singh, V. K., Kalsan, M., Kumar, N., Saini, A. & Chandra, R. Induced pluripotent stem cells: applications in regenerative medicine, disease modeling, and drug discovery. _Front Cell Dev Biol_ **3**, 2 (2015).

6.         New cells for old: The emerging potential of pluripotent stem cells in regenerative medicine.

7.         Nichols, J. & Smith, A. Naive and Primed Pluripotent States. _Cell Stem Cell_ **4**, 487–492 (2009).

8.         Kumari, D. States of Pluripotency: Naïve and Primed Pluripotent Stem Cells. in _Pluripotent Stem Cells - From the Bench to the Clinic_ (IntechOpen, 2016). doi:10.5772/63202.

9.         Smith, A. Formative pluripotency: the executive phase in a developmental continuum. _Development_ **144**, 365–373 (2017).

10.       Kinoshita, M. _et al._ Capture of Mouse and Human Stem Cells with Features of Formative Pluripotency. _Cell Stem Cell_ **28**, 453-471.e8 (2021).

11.       Shyh-Chang, N. & Li, L. Stabilizing Formative Pluripotent States with Germ Cell Competency. _Cell Stem Cell_ **28**, 361–363 (2021).

12.       Wang, X. _et al._ Formative pluripotent stem cells show features of epiblast cells poised for gastrulation. _Cell Res_ **31**, 526–541 (2021).

13.       Pera, M. F. & Rossant, J. The exploration of pluripotency space: Charting cell state transitions in peri-implantation development. _Cell Stem Cell_ **28**, 1896–1906 (2021).

14.       Ávila-González, D. _et al._ Unraveling the Spatiotemporal Human Pluripotency in Embryonic Development. _Front. Cell Dev. Biol._ **9**, (2021).

15.       Takahashi, S., Kobayashi, S. & Hiratani, I. Epigenetic differences between naïve and primed pluripotent stem cells. _Cell Mol Life Sci_ **75**, 1191–1203 (2018).

16.       Pelham-Webb, B., Murphy, D. & Apostolou, E. Dynamic 3D Chromatin Reorganization during Establishment and Maintenance of Pluripotency. _Stem Cell Reports_ **15**, 1176–1195 (2020).

17.       Platania, A. & Sexton, T. Chapter 4 - Chromatin architecture and topology in pluripotent stem cells. in _Stem Cell Epigenetics_ (eds. Meshorer, E. & Testa, G.) vol. 17 93–113 (Academic Press, 2020).

18.       Infrared and Raman Spectroscopic Imaging, 2nd, Completely Revised and Updated Edition | Wiley. _Wiley.com_ https://www.wiley.com/en-us/Infrared+and+Raman+Spectroscopic+Imaging%2C+2nd%2C+Completely+Revised+and+Updated+Edition-p-9783527336524.

19.       Imaging the invisible—Bioorthogonal Raman probes for imaging of cells and tissues - Azemtsop Matanfack - 2020 - Journal of Biophotonics - Wiley Online Library. https://onlinelibrary.wiley.com/doi/full/10.1002/jbio.202000129.

20.       Pérez-Jiménez, A. I., Lyu, D., Lu, Z., Liu, G. & Ren, B. Surface-enhanced Raman spectroscopy: benefits, trade-offs and future developments. _Chem. Sci._ **11**, 4563–4577 (2020).

21.       Gao, L. _et al._ Atomic Force Microscopy Based Tip-Enhanced Raman Spectroscopy in Biology. _Int J Mol Sci_ **19**, 1193 (2018).

22.       Chen, C., Hayazawa, N. & Kawata, S. A 1.7 nm resolution chemical analysis of carbon nanotubes by tip-enhanced Raman imaging in the ambient. _Nat Commun_ **5**, 3312 (2014).

23.       Treffer, R., Lin, X., Bailo, E., Deckert-Gaudig, T. & Deckert, V. Distinction of nucleobases - a tip-enhanced Raman approach. _Beilstein J Nanotechnol_ **2**, 628–637 (2011).

24.       Richter, M., Hedegaard, M., Deckert-Gaudig, T., Lampen, P. & Deckert, V. Laterally Resolved and Direct Spectroscopic Evidence of Nanometer-Sized Lipid and Protein Domains on a Single Cell. _Small_ **7**, 209–214 (2011).

25.       Zhao, Z. _et al._ Ultra-bright Raman dots for multiplexed optical imaging. _Nat Commun_ **12**, 1305 (2021).

26.       Chen, C. _et al._ Multiplexed live-cell profiling with Raman probes. _Nat Commun_ **12**, 3405 (2021).

27.       Optimized protocol for naive human pluripotent stem cell-derived trophoblast induction. https://star-protocols.cell.com/protocols/1135.

28.       Guo, G. _et al._ Epigenetic resetting of human pluripotency. _Development_ **144**, 2748–2763 (2017).

29.       Isolation of Nuclei for Single Cell RNA Sequencing & Tissues for Single Cell RNA Sequencing - Official 10x Genomics Support. _10x Genomics_ https://www.10xgenomics.com/support/single-cell-gene-expression/documentation/steps/sample-prep/isolation-of-nuclei-for-single-cell-rna-sequencing-and-tissues-for-single-cell-rna-sequencing.

30.       Branco, M. R. & Pombo, A. Intermingling of Chromosome Territories in Interphase Suggests Role in Translocations and Transcription-Dependent Associations. _PLOS Biology_ **4**, e138 (2006).

31.       Sanderson, M. J., Smith, I., Parker, I. & Bootman, M. D. Fluorescence Microscopy. _Cold Spring Harb Protoc_ **2014**, pdb.top071795 (2014).

32.       Galbraith, C. G. & Galbraith, J. A. Super-resolution microscopy at a glance. _J Cell Sci_ **124**, 1607–1611 (2011).

33.       Reinhardt, S. C. M. _et al._ Ångström-resolution fluorescence microscopy. _Nature_ **617**, 711–716 (2023).

34.       Koster, A. J. & Klumperman, J. Electron microscopy in cell biology: integrating structure and function. _Nat Rev Mol Cell Biol_ **Suppl**, SS6-10 (2003).

35.       Tillberg, P. W. Development of multiplexing strategies for electron and super-resolution optical microscopy/. (Massachusetts Institute of Technology, 2013).

36.       de Boer, P., Hoogenboom, J. P. & Giepmans, B. N. G. Correlated light and electron microscopy: ultrastructure lights up! _Nat Methods_ **12**, 503–513 (2015).

37.       Noble, D. It’s time to admit that genes are not the blueprint for life. _Nature_ **626**, 254–255 (2024).