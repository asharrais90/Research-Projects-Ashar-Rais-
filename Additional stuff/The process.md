- (thymus OR thymic OR "Thymus Gland"[Mesh] OR "Thymocytes"[Mesh] OR "thymic stroma") AND ("Biomechanical Phenomena"[Mesh] OR "Extracellular Matrix"[Mesh] OR "Mechanotransduction, Cellular"[Mesh] OR stiffness OR mechanobiology OR biomechanics OR mechanotransduction OR "mechanical properties" OR fibrosis OR involution OR "elastic modulus" OR "Young's modulus" OR rigidity OR viscoelastic OR mechanosens*) AND ("T-Lymphocytes"[Mesh] OR "Immune Tolerance"[Mesh] OR "Receptors, Antigen, T-Cell"[Mesh] OR thymocyte OR "T-cell selection" OR "positive selection" OR "negative selection" OR "agonist selection" OR tolerance OR "TCR signaling" OR "TCR signalling")
- ( TITLE-ABS-KEY ( thymus OR thymic OR "thymus gland" OR "thymic stroma" OR thymocytes OR "t-lymphocytes" OR "thymocyte development" OR "t-cell development" ) ) AND ( TITLE-ABS-KEY ( "biomechanical phenomena" OR "extracellular matrix" OR "mechanotransduction, cellular" OR stiffness OR mechanobiology OR biomechanics OR mechanotransduction OR "mechanical properties" OR fibrosis OR involution OR "elastic modulus" OR "young's modulus" OR rigidity OR viscoelastic OR mechanosens* ) ) AND ( TITLE-ABS-KEY ( "clonal selection, antigen-dependent" OR "immune tolerance" OR "receptors, antigen, t-cell" OR "t-cell selection" OR "positive selection" OR "negative selection" OR "agonist selection" OR "central tolerance" OR tolerance OR "tcr signalling" ) ) AND PUBYEAR > 2009 AND PUBYEAR < 2027 AND ( EXCLUDE ( DOCTYPE , "no" ) OR EXCLUDE ( DOCTYPE , "ed" ) OR EXCLUDE ( DOCTYPE , "sh" ) OR EXCLUDE ( DOCTYPE , "er" ) OR EXCLUDE ( DOCTYPE , "tb" ) OR EXCLUDE ( DOCTYPE , "le" ) ) AND ( LIMIT-TO ( LANGUAGE , "English" ) )
- - SM1: ("Thymus Gland"[Mesh] OR thymus[tiab] OR thymic[tiab] OR thymocyte*[tiab] OR "thymic microenvironment" [tiab] OR "thymic strom*"[tiab] OR "positive selection"[tiab] OR "negative selection"[tiab] OR "central tolerance"[tiab] OR "agonist selection" [tiab] OR "T cell receptor"[tiab] OR "T cell receptor signaling" [tiab] OR TCR[tiab] OR "TCR signaling" [tiab] OR "thymocyte differentiation" [tiab])
- SM2: ("Biomechanical Phenomena"[Mesh] OR "Extracellular Matrix"[Mesh] OR "Mechanotransduction, Cellular"[Mesh] OR viscoelasticity OR stiffness OR mechanobiology OR biomechanics OR mechanotransduction OR "mechanical properties" OR "matrix architecture" [tiab] OR "fiber orientation" [tiab] OR anisotropy [tiab] OR topography [tiab] OR collagen [tiab] OR Collagens [Mesh] OR laminin [tiab] OR Laminin [Mesh] OR fibronectin [tiab] OR Fibronectins [Mesh] OR "Integrins"[Mesh] OR integrin*[tiab] OR "Talin"[Mesh] OR talin[tiab] OR "Vinculin"[Mesh] OR vinculin[tiab] OR "Piezo1"[tiab]  OR YAP[tiab] OR TAZ[tiab] OR "Focal Adhesion Kinase 2"[Mesh] OR FAK[tiab] OR "rhoA GTP-Binding Protein"[Mesh]  OR "Actomyosin"[Mesh] OR actomyosin[tiab] OR "Cytoskeleton"[Mesh] OR cytoskeleton[tiab])
- SM3: ("Tissue Engineering"[Mesh] OR "Tissue Scaffolds"[Mesh] OR "Hydrogels"[Mesh] OR hydrogel*[tiab] OR scaffold*[tiab] OR "thymus organoid*"[tiab] OR biomimetic[tiab] OR "engineered ECM"[tiab] OR "artificial thymus"[tiab] OR "3D culture"[tiab] OR decellular*[tiab] OR "Engineered thymic niches" [tiab])
- SM4: (involution[tiab] OR aging[tiab] OR ageing[tiab] OR fibrosis[tiab] OR "fibrotic remodeling" [tiab] OR "ECM remodeling" [tiab] OR "extracellular remodeling" [tiab] OR senescence[tiab])
- Included: Books and Documents, Clinical Study, Comparative Study, Dataset, Evidence Synthesis, Meta-Analysis, Observational Study, Preprint, Review, Scoping Review, Systematic Review, Validation Study.
Scopus: 3612 (Search 1) results + 22 pre-prints
- SM1: TITLE-ABS-KEY( "Thymus Gland" OR thymus OR thymic OR thymocyte* OR "thymic microenvironment" OR "thymic stroma" OR "thymic stromal" OR "positive selection" OR "negative selection" OR "central tolerance" OR "agonist selection" OR "T cell receptor" OR "T cell receptor signaling" OR TCR OR "TCR signaling" OR "thymocyte differentiation")
- SM2: TITLE-ABS-KEY( "Biomechanical Phenomena" OR "Extracellular Matrix" OR "Mechanotransduction, Cellular" OR viscoelasticity OR stiffness OR mechanobiology OR biomechanics OR mechanotransduction OR "mechanical properties" OR "matrix architecture" OR "fiber orientation" OR anisotropy OR topography OR Collagen OR Laminin OR Fibronectin OR integrin* OR Talin OR VinculiN OR Piezo1 OR YAP OR TAZ OR "Focal Adhesion Kinase 2" OR FAK OR "rhoA GTP-Binding Protein" OR Actomyosin OR Cytoskeleton)
- SM3: TITLE-ABS-KEY( "Tissue Engineering" OR "Tissue Scaffolds" OR "Hydrogels" OR hydrogel* OR scaffold* OR "thymus organoid" OR "thymus organoids" OR biomimetic OR "engineered ECM" OR "artificial thymus" OR "3D culture" OR decellular* OR "Engineered thymic niches")
- SM4: TITLE-ABS-KEY( involution OR aging OR ageing OR fibrosis OR "fibrotic remodeling" OR "ECM remodeling" OR "extracellular remodeling" OR senescence)
- - After deduplication: 977 items;
- After title screening: 92 items (9.4%);
- After abstract screening: 35 items (38%);
- After full-text screening: 16 items (~50%).



Questions for Claude:
1.     def _S(self) -> float:

        x, y = self.pos

        return float(self.model.stiffness[x, y]). How does this work and would there be a difference if we typed self.pos = x,y
2. def _E(self) -> int:

        """E(t)=1 if any Moore-neighbour is a stromal cell."""

        for nb in self.model.grid.iter_neighborhood(

                self.pos, moore=True, include_center=False):

            if nb in self.model.stromal_pos:

                return 1

        return 0
      Why is there no 'else' here?
  3. Why is p = self.model.p used instead of p = self.p. Since a separate parameter class has been defined, wouldn't it be fine even if each thymocyte accesses p directly instead of attaching p to the model and then having the agents access the values? Same question for rng = self.model.rng (Also, is this equivalent to the critique you answered a little while back?)?
  4.  Can you explain this code line-by-line please:         rising  = S**p.n_low  / (p.K_low**p.n_low  + S**p.n_low)

        falling = p.K_high**p.n_high  / (p.k_high**p.n_high + S**p.n_high)

        if rng.random() > rising * falling:

            return

  

        candidates = list(self.model.grid.iter_neighborhood(

            self.pos, moore=True, include_center=False))

        if not candidates:

            return

  

        if rng.random() < p.p_bias:

            biased = [c for c in candidates if c[1] > self.pos[1]]

            if biased:

                candidates = biased

  

        chosen = candidates[rng.integers(len(candidates))]

        if self.model.grid.is_cell_empty(chosen):

            self.model.grid.move_agent(self, chosen)
        Also, the value of rising and falling are bound between 0 and 1 right? Also, where is v_max?
        Why is list explicitly mentioned, Also, correct me if I am wrong: the candidates variable (which contains the list of grid points around each agent) is formed regardless of p_bias. And another variable biased is formed which checks the y coordinates of the grid points (why is 1 used instead of 0,1) and if the y-coordinate of any point is greater than the y-coordinate of the agent, the candidate variable is updated to contain only the biased points (so those with > y-coordinates). Now you create another variable 'chosen' which takes a random value from the candidates list. Finally the method 'self.model.grid.is_cell_empty' is called and if this is true, then the move_agent method is called with the arguments self and chosen (so, move this particular agent to the chosen grid point). Also, why is self passed as an argument, when it is already mentioned in the calling of the method? 


        5. self.C = float(np.clip(self.C + self._dC_dt(S)    * dt, 0.0, 1.0))

        self.M = float(np.clip(self.M + self._dM_dt(S)    * dt, 0.0, 1.0))

        self.X = max(0.0,      self.X + self._dX_dt(E, S) * dt)_
       Why doesn't self.X has the float conversion? Or maybe the right question is that why do C and M have this conversion?
6.  self.rng = np.random.default_rng(self.p.seed). What is the default rng?
7. Outputs to show the average velocity and migrational trajectories (in a stiffness map) of thymocytes in each of the outcome groups. 
8. Running the simulation until the last agent has either negatively selected or completed its residence time.
9. Saturation of C, M, k_scan and alpha. And if we see a spike in X at a later time what kind of useful information might be extrapolated from this and what additional functionalities should be added to make this graph more useful?
10. what is "Window width: {np.log10(tn/tp):.2f} log₁₀ decades"


Assalamu 'alaikum,

Thanks for going through this in detail, it kind of exceeded my expectations.

And I should have sent the updated rules with the code, but my main intent was to decrease the chance of having a major computational hole in the simulation (Also, I was rushing things a bit because I was/still am behind my schedule). 

Anyways, I am attaching the improved version of the code if you want to take a look, test a few things and give your opinion about its quality, unchecked mistakes, etc. Also, just to let you know my plan with producing results from this simple code is to generate some falsifiable hypotheses that might (along with the actual work) interest a lab-head who is working on similar topics (I have already found many people whom I can contact). One of the main things about this super-specialized topic is that there are many, many things that haven't been explored. Most importantly for this simulation, a consensus hasn't been established yet about the stiffness range of healthy thymus stroma (many studies producing different results like 1-3 and 6-12 kPa depending on the samples and techniques used), whether the cortex is stiffer or the medulla and how these biomechanical factors influence thymic developmental steps. Consequently, I plan on running this model with different ranges (1-3, 3-6, 6-12, 15-20, etc.) each with 3 scenarios (the uniform, cortex stiffer...) to try and get some useful results (for eg, if the cortex_stiffer scenario produces results matching real-life experimental data on relative fractions of selection outcomes for most stiffness ranges while the medulla_stiffere scenario doesn't produce consistent results than this might be a possibly useful datapoint in reaching consensus about this particular question).  
Also to answer some of your thoughtful comments:  
1 (i). At first I included the mechanotransduction module that connects directly to the fate outcomes (hence the "mechanically rescued" term that you might have seen in the v2 that should have been removed if I wasn't rushing things so much) because pathways such as YAP/TAZ and MRTF/SRF have been observed to influence cell survival across different cell types. However, the main issue with this was that I did not have any experiment-based way to connect their contribution to the experimentally backed TCR-signalling thresholds and doing this arbitrarily seemed imprudent. But during my parallel work on the scoping review, I read a few papers which mentione the role of mechanotransductive receptors like Piezo-1 in influencing TCR signalling, so I decided to take a different route for bringing these mechanotransductive processes into the model;

(ii) Yeah, that was decided after I sent the first rule draft to you (again, due to reading some papers for the review);

2. As this an estimated parameter (based on real-life data) I am planning to include it in the sensitivity analysis, although this won't have as much as priority as given to it when compared to other phenomenological parameters like k_scan,base, alpha_o, etc.;

3. Thank you point that out I have made the relevant change:         for nb in self.model.grid.iter_neighborhood(

          self.pos, moore=True, include_center=True):

6.