export const cellularComponents = {
    id: 'cellular-components',
    name: 'Cellular Components',
    description: 'Fundamental building blocks and molecular machinery that constitute the structural and functional units of single cells',
    
    categories: [
        {
            id: 'membranes',
            name: 'Membranes & Barriers',
            icon: 'fas fa-shield-alt',
            color: 'primary',
            items: [
                {
                    id: 'plasma-membrane',
                    name: 'Plasma Membrane',
                    description: 'The selectively permeable lipid bilayer that defines the cell boundary and mediates interactions with the environment',
                    components: [
                        'Phospholipid bilayer (40-50% of mass)',
                        'Integral membrane proteins (25-30%)',
                        'Peripheral membrane proteins (5%)',
                        'Cholesterol (20-25% of lipids)',
                        'Glycocalyx (carbohydrate coat)'
                    ],
                    composition: {
                        'Phospholipids': 'PC (45%), PE (25%), PS (10%), PI (10%), others',
                        'Sphingolipids': 'Sphingomyelin, glycosphingolipids',
                        'Cholesterol': '20-50% of membrane lipids',
                        'Proteins': '~10,000 different types per cell'
                    },
                    asymmetry: {
                        'Outer leaflet': 'PC, sphingomyelin, glycolipids',
                        'Inner leaflet': 'PS, PE, PI, PIP2',
                        'Maintenance': 'ATP-dependent flippases'
                    },
                    techniques: [
                        'Patch-clamp electrophysiology',
                        'FRAP (Fluorescence Recovery After Photobleaching)',
                        'Single particle tracking',
                        'Membrane proteomics',
                        'Lipidomics analysis',
                        'Super-resolution microscopy (STORM, PALM)'
                    ],
                    databases: [
                        { name: 'TCDB', url: 'https://www.tcdb.org/', description: 'Transporter Classification Database' },
                        { name: 'OPM', url: 'https://opm.phar.umich.edu/', description: 'Orientations of Proteins in Membranes' },
                        { name: 'MemProtMD', url: 'http://memprotmd.bioch.ox.ac.uk/', description: 'Membrane protein structures' },
                        { name: 'PDBTM', url: 'https://pdbtm.unitmp.org/', description: 'Protein Data Bank of Transmembrane Proteins' }
                    ],
                    keyFacts: [
                        'Thickness: 7-10 nm',
                        'Surface area: ~700 µm² (typical mammalian cell)',
                        'Lipid molecules: ~5×10^8 per cell',
                        'Lateral diffusion: 10^-8 cm²/s for lipids',
                        'Protein density: ~20,000 proteins/µm²'
                    ],
                    functions: [
                        'Selective permeability barrier',
                        'Signal transduction platform',
                        'Cell adhesion and recognition',
                        'Energy transduction (gradients)',
                        'Endocytosis and exocytosis'
                    ]
                },
                {
                    id: 'nuclear-envelope',
                    name: 'Nuclear Envelope',
                    description: 'Double membrane system that compartmentalizes the genome and regulates nucleocytoplasmic transport',
                    structure: {
                        'Outer membrane': 'Continuous with rough ER, ribosomes attached',
                        'Inner membrane': 'Unique proteins, lamina attachment',
                        'Perinuclear space': '20-40 nm width',
                        'Nuclear pores': '~2000-4000 per nucleus',
                        'Nuclear lamina': 'Lamin A/C, B1, B2 meshwork'
                    },
                    nuclearPoreComplex: {
                        'Size': '~125 MDa, 120 nm diameter',
                        'Composition': '~30 different nucleoporins',
                        'Transport rate': '~1000 translocations/pore/second',
                        'Size exclusion': '~40 kDa passive, >40 kDa active',
                        'FG-Nups': 'Intrinsically disordered barrier'
                    },
                    techniques: [
                        'Super-resolution microscopy',
                        'Nuclear import/export assays',
                        'FRAP for transport dynamics',
                        'Electron microscopy',
                        'DamID for lamina interactions'
                    ],
                    // databases: [
                    //     { name: 'NPC Database', url: 'http://www.npcdb.org/', description: 'Nuclear Pore Complex components' },
                    //     { name: 'NuPoP', url: 'https://www.predictNLS.org/', description: 'Nuclear localization signals' }
                    // ],
                    diseases: [
                        'Laminopathies (LMNA mutations)',
                        'Hutchinson-Gilford progeria',
                        'Emery-Dreifuss muscular dystrophy',
                        'Nuclear envelope spectrin repeat proteins (Nesprins) disorders'
                    ]
                },
                {
                    id: 'er-membrane',
                    name: 'Endoplasmic Reticulum Membrane',
                    description: 'Extensive membrane network for protein synthesis, folding, and lipid biosynthesis',
                    organization: {
                        'Rough ER': 'Ribosome-studded, protein synthesis',
                        'Smooth ER': 'Lipid synthesis, Ca²⁺ storage',
                        'Transitional ER': 'COPII vesicle formation',
                        'MAMs': 'Mitochondria-associated membranes'
                    },
                    proteins: [
                        'Sec61 translocon complex',
                        'OST (oligosaccharyltransferase)',
                        'Protein disulfide isomerase (PDI)',
                        'BiP/GRP78 chaperone',
                        'Calnexin/Calreticulin'
                    ],
                    lipidSynthesis: [
                        'Phospholipid biosynthesis',
                        'Cholesterol synthesis enzymes',
                        'Ceramide synthesis',
                        'GPI anchor assembly'
                    ],
                    databases: [
                        { name: 'ERnet', url: 'http://www.secretomes.org/', description: 'ER and secretory pathway' }
                    ]
                },
                {
                    id: 'mitochondrial-membranes',
                    name: 'Mitochondrial Membranes',
                    description: 'Dual membrane system enabling ATP production and metabolic compartmentalization',
                    structure: {
                        'Outer membrane': 'Permeable to <5 kDa, VDAC channels',
                        'Inner membrane': 'Highly folded cristae, ~75% protein',
                        'Intermembrane space': 'Cytochrome c, creatine kinase',
                        'Matrix': 'mtDNA, ribosomes, metabolic enzymes'
                    },
                    complexes: [
                        'Complex I: NADH dehydrogenase (45 subunits)',
                        'Complex II: Succinate dehydrogenase',
                        'Complex III: Cytochrome bc1',
                        'Complex IV: Cytochrome c oxidase',
                        'Complex V: ATP synthase'
                    ],
                    dynamics: {
                        'Fusion': 'Mfn1/2 (outer), OPA1 (inner)',
                        'Fission': 'Drp1, Fis1, MFF',
                        'Cristae remodeling': 'OPA1, MICOS complex',
                        'Mitophagy': 'PINK1/Parkin pathway'
                    },
                    databases: [
                        { name: 'MitoCarta', url: 'https://www.broadinstitute.org/mitocarta/mitocarta30-inventory-mammalian-mitochondrial-proteins-and-pathways', description: 'Mitochondrial protein inventory' },
                        { name: 'MitoMiner', url: 'http://mitominer.mrc-mbu.cam.ac.uk/', description: 'Mitochondrial proteomics database' }
                    ]
                }
            ]
        },
        {
            id: 'organelles',
            name: 'Organelles',
            icon: 'fas fa-shapes',
            color: 'secondary',
            items: [
                {
                    id: 'nucleus',
                    name: 'Nucleus',
                    description: 'The control center containing genetic material and coordinating cellular activities including growth, metabolism, and reproduction',
                    substructures: [
                        'Nucleolus - ribosome biogenesis, ~1-5 per nucleus',
                        'Nuclear speckles - splicing factor storage',
                        'Cajal bodies - snRNP assembly',
                        'PML bodies - transcriptional regulation',
                        'Chromatin territories - chromosome domains',
                        'Nuclear matrix - structural scaffold'
                    ],
                    chromatin: {
                        'Euchromatin': 'Active, accessible, gene-rich',
                        'Heterochromatin': 'Condensed, silenced, repetitive',
                        'Nucleosome': '147 bp DNA + histone octamer',
                        'Chromosome territories': 'Non-random 3D organization',
                        'TADs': 'Topologically associating domains'
                    },
                    singleCellMethods: [
                        'scATAC-seq - chromatin accessibility',
                        'scHi-C - 3D genome organization',
                        'DNA FISH - specific loci visualization',
                        'scDam-ID - lamina interactions',
                        'CUT&Tag - histone modifications'
                    ],
                    databases: [
                        { name: 'Nuclear Protein Database', url: 'https://www.rostlab.org/services/npd/', description: 'Nuclear proteome' },
                        { name: 'NucleomeDB', url: 'http://nucleomedb.org/', description: 'Nuclear proteome database' },
                        { name: '4DN Portal', url: 'https://data.4dnucleome.org/', description: '4D nucleome organization' },
                        { name: 'ChromatinDB', url: 'http://www.bio.ifi.lmu.de/ChromatinDB/', description: 'Chromatin-associated proteins' }
                    ],
                    keyMetrics: {
                        'Volume': '~520 µm³ (HeLa cells)',
                        'DNA content': '6.4 Gb diploid human genome',
                        'Genes': '~20,000 protein-coding',
                        'Nuclear pores': '2000-4000 per nucleus',
                        'Transcription factories': '~100-300 active sites'
                    },
                    nuclearBodies: [
                        'Nucleolus (1-5 per cell, ~1-3 µm)',
                        'Speckles (20-50 per cell)',
                        'Cajal bodies (1-10 per cell)',
                        'PML bodies (10-30 per cell)',
                        'Paraspeckles (2-20 per cell)'
                    ]
                },
                {
                    id: 'mitochondria',
                    name: 'Mitochondria',
                    description: 'Double-membrane organelles that generate ATP through oxidative phosphorylation and regulate cellular metabolism',
                    features: [
                        'Double membrane structure',
                        'Own circular DNA genome (16.5 kb, 37 genes)',
                        'Semi-autonomous replication',
                        'Dynamic network through fusion/fission',
                        'Maternal inheritance pattern'
                    ],
                    functions: {
                        'Energy': 'ATP production via OXPHOS',
                        'Metabolism': 'TCA cycle, β-oxidation, amino acids',
                        'Signaling': 'ROS production, Ca²⁺ buffering',
                        'Biosynthesis': 'Heme, iron-sulfur clusters, steroids',
                        'Cell death': 'Apoptosis regulation (cytochrome c)'
                    },
                    singleCellAspects: [
                        'Heteroplasmy detection in scDNA-seq',
                        'Mitochondrial gene expression in scRNA-seq',
                        'Quality metrics (% mitochondrial reads)',
                        'Metabolic state inference',
                        'Lineage tracing via mutations'
                    ],
                    dynamics: {
                        'Number': '100-2000 per cell type-dependent',
                        'Volume': '~20% of cell volume',
                        'Turnover': 'Half-life ~10-25 days',
                        'Movement': 'Along microtubules, ~0.5 µm/s'
                    },
                    databases: [
                        { name: 'MitoMiner', url: 'http://mitominer.mrc-mbu.cam.ac.uk/', description: 'Integrated mitochondrial proteomics' },
                        { name: 'MitoCarta3.0', url: 'https://www.broadinstitute.org/mitocarta/', description: 'Mitochondrial protein inventory' },
                        { name: 'MITOMAP', url: 'https://www.mitomap.org/', description: 'Human mitochondrial genome database' },
                        { name: 'MitoAge', url: 'http://www.mitoage.info/', description: 'Mitochondria and aging' }
                    ]
                },
                {
                    id: 'endoplasmic-reticulum',
                    name: 'Endoplasmic Reticulum',
                    description: 'Extensive membrane network for protein synthesis, folding, modification, and lipid biosynthesis',
                    types: {
                        'Rough ER': {
                            function: 'Protein synthesis and folding',
                            features: 'Ribosome-studded surface',
                            products: 'Secretory and membrane proteins'
                        },
                        'Smooth ER': {
                            function: 'Lipid synthesis, detoxification',
                            features: 'No ribosomes',
                            products: 'Lipids, steroids, Ca²⁺ storage'
                        }
                    },
                    functions: [
                        'Co-translational protein insertion',
                        'N-glycosylation initiation',
                        'Disulfide bond formation',
                        'Quality control (ERAD)',
                        'Lipid biosynthesis',
                        'Calcium homeostasis'
                    ],
                    markers: [
                        'BiP/GRP78 - ER chaperone',
                        'Calnexin - glycoprotein folding',
                        'PDI - disulfide isomerase',
                        'Sec61 - translocon',
                        'KDEL receptor - retrieval'
                    ],
                    stress: {
                        'UPR': 'Unfolded protein response',
                        'Sensors': 'IRE1, PERK, ATF6',
                        'Outcomes': 'Adaptation, autophagy, apoptosis'
                    },
                    databases: [
                        { name: 'Human Protein Atlas - ER', url: 'https://www.proteinatlas.org/subcellular/endoplasmic+reticulum', description: 'ER protein localization' }
                    ]
                },
                {
                    id: 'golgi',
                    name: 'Golgi Apparatus',
                    description: 'Protein and lipid processing, modification, and sorting center for the secretory pathway',
                    structure: [
                        'Cis-Golgi network (CGN) - receiving',
                        'Cis-cisternae - early processing',
                        'Medial-cisternae - core modifications',
                        'Trans-cisternae - late processing',
                        'Trans-Golgi network (TGN) - sorting'
                    ],
                    modifications: [
                        'N-glycan processing',
                        'O-glycosylation',
                        'Proteolytic cleavage',
                        'Sulfation',
                        'Phosphorylation of lysosomal proteins'
                    ],
                    transport: {
                        'Anterograde': 'ER to Golgi via COPII',
                        'Retrograde': 'Golgi to ER via COPI',
                        'Intra-Golgi': 'Cisternal maturation model',
                        'Exit': 'Clathrin-coated vesicles from TGN'
                    },
                    markers: [
                        'GM130 - cis-Golgi matrix',
                        'Giantin - Golgi membrane',
                        'TGN46 - trans-Golgi network',
                        'Golgin-97 - TGN',
                        'GRASP65/55 - stacking'
                    ],
                    databases: [
                        { name: 'GlycoStore', url: 'https://glycostore.org/', description: 'Glycan structures database' }
                    ]
                },
                {
                    id: 'lysosomes',
                    name: 'Lysosomes',
                    description: 'Acidic organelles containing hydrolytic enzymes for cellular digestion and recycling',
                    properties: {
                        'pH': '4.5-5.0 maintained by V-ATPase',
                        'Enzymes': '~60 different hydrolases',
                        'Membrane': 'Highly glycosylated protection',
                        'Number': '50-1000 per cell',
                        'Size': '0.1-1.2 µm diameter'
                    },
                    enzymes: [
                        'Proteases (cathepsins)',
                        'Nucleases',
                        'Lipases',
                        'Glycosidases',
                        'Phosphatases',
                        'Sulfatases'
                    ],
                    functions: [
                        'Autophagy',
                        'Endocytosis degradation',
                        'Pathogen destruction',
                        'Membrane repair',
                        'Cell signaling (mTORC1)',
                        'Cholesterol homeostasis'
                    ],
                    diseases: [
                        'Lysosomal storage diseases (>50 types)',
                        'Gaucher disease',
                        'Pompe disease',
                        'Niemann-Pick disease'
                    ],
                    markers: [
                        'LAMP1/2 - membrane proteins',
                        'Cathepsin D - protease',
                        'V-ATPase - proton pump'
                    ]
                },
                {
                    id: 'peroxisomes',
                    name: 'Peroxisomes',
                    description: 'Single-membrane organelles involved in lipid metabolism and reactive oxygen species regulation',
                    functions: [
                        'β-oxidation of very long chain fatty acids',
                        'Plasmalogen synthesis',
                        'Bile acid synthesis',
                        'H₂O₂ metabolism (catalase)',
                        'Glyoxylate detoxification'
                    ],
                    biogenesis: {
                        'De novo': 'From ER-derived vesicles',
                        'Fission': 'Division of existing peroxisomes',
                        'Import': 'PTS1/PTS2 targeting signals'
                    },
                    markers: [
                        'PMP70 - membrane protein',
                        'Catalase - H₂O₂ breakdown',
                        'PEX proteins - biogenesis'
                    ],
                    diseases: [
                        'Zellweger syndrome',
                        'Adrenoleukodystrophy',
                        'Refsum disease'
                    ]
                }
            ]
        },
        {
            id: 'molecular-machinery',
            name: 'Molecular Machinery',
            icon: 'fas fa-cogs',
            color: 'success',
            items: [
                {
                    id: 'ribosomes',
                    name: 'Ribosomes',
                    description: 'Molecular machines that translate mRNA into proteins, consisting of ribosomal RNA and proteins',
                    composition: {
                        'Prokaryotic': '70S (50S + 30S subunits)',
                        'Eukaryotic': '80S (60S + 40S subunits)',
                        'rRNA': '~60% of mass, catalytic',
                        'Proteins': '~40% of mass, structural',
                        'Size': '~25-30 nm diameter'
                    },
                    types: [
                        'Free ribosomes - cytosolic proteins',
                        'ER-bound ribosomes - secretory/membrane proteins',
                        'Mitochondrial ribosomes - 55S mitoribosomes',
                        'Plastid ribosomes - 70S in chloroplasts'
                    ],
                    dynamics: {
                        'Translation rate': '10-20 amino acids/second',
                        'Number per cell': '~10 million in mammalian cells',
                        'Polysome formation': 'Multiple ribosomes per mRNA',
                        'Turnover': 'Half-life ~5 days'
                    },
                    singleCellRelevance: [
                        'Translation efficiency varies between cells',
                        'Ribosome profiling (Ribo-seq)',
                        'Ribosome heterogeneity',
                        'Specialized ribosomes concept'
                    ],
                    databases: [
                        { name: 'RCSB PDB', url: 'https://www.rcsb.org/', description: 'Ribosome structures' },
                        { name: 'Ribosome Database', url: 'https://ribosome.med.nyu.edu/', description: 'Ribosomal RNA sequences' }
                    ]
                },
                {
                    id: 'proteasome',
                    name: 'Proteasome',
                    description: 'Multi-subunit protease complex responsible for regulated protein degradation',
                    structure: {
                        '26S proteasome': '20S core + 19S regulatory (2.5 MDa)',
                        '20S core': '4 stacked rings (α₇β₇β₇α₇)',
                        '19S regulatory': 'Recognition and unfolding',
                        'Immunoproteasome': 'Specialized for antigen processing',
                        'Thymoproteasome': 'T-cell selection'
                    },
                    mechanism: [
                        'Ubiquitin recognition',
                        'Deubiquitination',
                        'Substrate unfolding',
                        'Translocation into core',
                        'Proteolytic cleavage',
                        'Peptide release'
                    ],
                    regulation: {
                        'Ubiquitin code': 'K48, K63, K11 linkages',
                        'E3 ligases': '>600 in humans',
                        'DUBs': '~100 deubiquitinases',
                        'Proteasome assembly': 'Dedicated chaperones'
                    },
                    inhibitors: [
                        'Bortezomib (Velcade) - myeloma treatment',
                        'Carfilzomib - second generation',
                        'MG132 - research tool',
                        'Epoxomicin - natural product'
                    ],
                    databases: [
                        { name: 'MEROPS', url: 'https://www.ebi.ac.uk/merops/', description: 'Protease database' }
                    ]
                },
                {
                    id: 'cytoskeleton',
                    name: 'Cytoskeleton',
                    description: 'Dynamic network of protein filaments providing structure, transport tracks, and force generation',
                    components: [
                        {
                            name: 'Microfilaments (Actin)',
                            protein: 'G-actin → F-actin',
                            diameter: '7 nm',
                            dynamics: 'Treadmilling, branching (Arp2/3)',
                            functions: 'Cell shape, migration, division',
                            motors: 'Myosins (>40 types)'
                        },
                        {
                            name: 'Intermediate filaments',
                            proteins: 'Keratins, vimentin, lamins, neurofilaments',
                            diameter: '10 nm',
                            dynamics: 'Less dynamic, mechanical support',
                            functions: 'Mechanical strength, nuclear structure'
                        },
                        {
                            name: 'Microtubules',
                            protein: 'α/β-tubulin heterodimers',
                            diameter: '25 nm',
                            dynamics: 'Dynamic instability, catastrophe/rescue',
                            functions: 'Transport, mitosis, cilia',
                            motors: 'Kinesins (45 types), dyneins'
                        }
                    ],
                    regulation: [
                        'Small GTPases (RhoA, Rac1, Cdc42)',
                        'Actin-binding proteins (>100 types)',
                        'MAPs (microtubule-associated proteins)',
                        'Post-translational modifications'
                    ],
                    imaging: [
                        'Phalloidin - F-actin staining',
                        'Anti-tubulin antibodies',
                        'Live-cell probes (LifeAct, SiR-actin)',
                        'Super-resolution microscopy'
                    ],
                    databases: [
                        { name: 'CytoSkeleton Database', url: 'http://www.cytoskeleton.info/', description: 'Cytoskeletal proteins' }
                    ]
                },
                {
                    id: 'motor-proteins',
                    name: 'Molecular Motors',
                    description: 'ATP-powered proteins that generate force and movement along cytoskeletal tracks',
                    types: {
                        'Myosins': {
                            track: 'Actin filaments',
                            direction: 'Mostly plus-end (barbed)',
                            families: '>40 classes',
                            functions: 'Muscle contraction, transport, tension'
                        },
                        'Kinesins': {
                            track: 'Microtubules',
                            direction: 'Mostly plus-end (cell periphery)',
                            families: '45 in humans',
                            functions: 'Organelle transport, mitosis'
                        },
                        'Dyneins': {
                            track: 'Microtubules',
                            direction: 'Minus-end (cell center)',
                            types: 'Cytoplasmic and axonemal',
                            functions: 'Retrograde transport, cilia beating'
                        }
                    },
                    mechanics: {
                        'Step size': '8-36 nm',
                        'Force': '1-7 pN',
                        'Speed': '0.1-60 µm/s',
                        'Processivity': '1-100+ steps',
                        'Efficiency': '~50% ATP to work'
                    },
                    cargo: [
                        'Vesicles and organelles',
                        'mRNA-protein complexes',
                        'Chromosomes during mitosis',
                        'Signaling complexes'
                    ]
                },
                {
                    id: 'chaperones',
                    name: 'Molecular Chaperones',
                    description: 'Proteins that assist in proper folding and prevent aggregation of other proteins',
                    families: {
                        'Hsp70': {
                            function: 'Co-translational folding',
                            cofactors: 'J-proteins (Hsp40), NEFs',
                            ATP: 'Required for substrate cycling'
                        },
                        'Hsp90': {
                            function: 'Client protein maturation',
                            clients: 'Kinases, receptors, transcription factors',
                            cochaperones: '>20 regulatory proteins'
                        },
                        'Hsp60/GroEL': {
                            function: 'Folding chamber',
                            structure: 'Double-ring barrel',
                            location: 'Mitochondria, cytosol'
                        },
                        'Small HSPs': {
                            function: 'Prevent aggregation',
                            ATP: 'Independent',
                            examples: 'αB-crystallin, Hsp27'
                        }
                    },
                    stress: {
                        'Heat shock response': 'HSF1 activation',
                        'ER stress': 'BiP/GRP78, calnexin',
                        'Oxidative stress': 'Increased expression',
                        'Proteotoxic stress': 'Inclusion body formation'
                    },
                    diseases: [
                        'Protein misfolding diseases',
                        'Neurodegeneration',
                        'Cancer (Hsp90 clients)',
                        'Chaperonopathies'
                    ]
                }
            ]
        },
        {
            id: 'molecular-complexes',
            name: 'Molecular Complexes',
            icon: 'fas fa-project-diagram',
            color: 'warning',
            items: [
                {
                    id: 'transcription-machinery',
                    name: 'Transcription Machinery',
                    description: 'Multi-subunit complexes that synthesize RNA from DNA templates',
                    components: {
                        'RNA Pol II': {
                            subunits: '12 subunits, ~550 kDa',
                            CTD: 'C-terminal domain with 52 repeats',
                            function: 'mRNA, lncRNA, miRNA transcription'
                        },
                        'General TFs': {
                            factors: 'TFIIA, B, D, E, F, H',
                            function: 'Basal transcription initiation',
                            assembly: 'Pre-initiation complex formation'
                        },
                        'Mediator': {
                            subunits: '~30 subunits, 1.4 MDa',
                            modules: 'Head, middle, tail, kinase',
                            function: 'Bridge enhancers to promoters'
                        }
                    },
                    regulation: [
                        'Enhancers and promoters',
                        'Chromatin remodeling complexes',
                        'Histone modifications',
                        'Transcription factors (~1600 in humans)',
                        'Non-coding RNA regulation'
                    ],
                    elongation: {
                        'Speed': '1-4 kb/min',
                        'Pausing': 'Promoter-proximal, gene body',
                        'Factors': 'P-TEFb, DSIF, NELF',
                        'Coupling': 'Co-transcriptional processing'
                    },
                    singleCellMethods: [
                        'scRNA-seq - output measurement',
                        'scATAC-seq - accessibility',
                        'ChIP-seq - factor binding',
                        'PRO-seq - active transcription'
                    ],
                    databases: [
                        { name: 'JASPAR', url: 'https://jaspar.genereg.net/', description: 'Transcription factor binding profiles' },
                        { name: 'TRANSFAC', url: 'http://genexplain.com/transfac/', description: 'TF database' }
                    ]
                },
                {
                    id: 'spliceosome',
                    name: 'Spliceosome',
                    description: 'Dynamic ribonucleoprotein complex that removes introns from pre-mRNA',
                    components: [
                        'U1 snRNP - 5\' splice site recognition',
                        'U2 snRNP - branch point binding',
                        'U4/U6•U5 tri-snRNP - catalytic core',
                        '>150 proteins involved',
                        '~45 core spliceosomal proteins'
                    ],
                    assembly: {
                        'E complex': 'Early, U1 binding',
                        'A complex': 'U2 addition',
                        'B complex': 'Tri-snRNP joining',
                        'C complex': 'Catalytically active',
                        'P complex': 'Post-catalytic'
                    },
                    alternativeSplicing: {
                        'Prevalence': '>95% of human genes',
                        'Types': 'Exon skipping, alternative 3\'/5\' sites',
                        'Regulation': 'SR proteins, hnRNPs',
                        'Tissue-specific': 'Brain highest diversity'
                    },
                    diseases: [
                        'Spinal muscular atrophy (SMN)',
                        'Retinitis pigmentosa',
                        'Myelodysplastic syndromes',
                        'Many cancers (SF3B1 mutations)'
                    ],
                    databases: [
                        { name: 'SpliceDB', url: 'https://www.splice-db.org/', description: 'Alternative splicing database' }
                    ]
                },
                {
                    id: 'phase-separated-condensates',
                    name: 'Biomolecular Condensates',
                    description: 'Membrane-less organelles formed by liquid-liquid phase separation',
                    examples: [
                        'Nucleolus - ribosome biogenesis',
                        'Nuclear speckles - splicing factors',
                        'Stress granules - mRNA storage',
                        'P-bodies - mRNA decay',
                        'Cajal bodies - snRNP assembly',
                        'Germ granules - germline determinants'
                    ],
                    properties: {
                        'Formation': 'Weak multivalent interactions',
                        'Dynamics': 'Liquid-like, rapid exchange',
                        'Components': 'IDRs, RNA, scaffolds',
                        'Regulation': 'PTMs, RNA, ATP',
                        'Size': '0.1-10 µm typically'
                    },
                    principles: [
                        'Intrinsically disordered regions (IDRs)',
                        'Multivalent interactions',
                        'Concentration-dependent assembly',
                        'ATP-dependent regulation',
                        'RNA as scaffold or client'
                    ],
                    detection: [
                        'FRAP - liquid-like behavior',
                        '1,6-hexanediol - disruption test',
                        'Optogenetic assembly (optoDroplets)',
                        'Time-lapse microscopy',
                        'In vitro reconstitution'
                    ],
                    diseases: [
                        'ALS/FTD - FUS, TDP-43 aggregation',
                        'Cancer - oncogenic fusions',
                        'Repeat expansion diseases',
                        'Ribosomopathies'
                    ],
                    databases: [
                        { name: 'PhaSepDB', url: 'http://db.phasep.pro/', description: 'Phase separation database' },
                        { name: 'DrLLPS', url: 'http://llps.biocuckoo.cn/', description: 'Liquid-liquid phase separation' }
                    ]
                },
                {
                    id: 'dna-repair-complexes',
                    name: 'DNA Repair Complexes',
                    description: 'Multi-protein assemblies that detect and repair DNA damage',
                    pathways: {
                        'Base excision repair': {
                            proteins: 'DNA glycosylases, APE1, Pol β',
                            lesions: 'Oxidized bases, small adducts'
                        },
                        'Nucleotide excision repair': {
                            proteins: 'XPC, TFIIH, XPG, XPF',
                            lesions: 'UV damage, bulky adducts'
                        },
                        'Mismatch repair': {
                            proteins: 'MSH2/6, MLH1, PMS2',
                            lesions: 'Replication errors'
                        },
                        'Double-strand break repair': {
                            NHEJ: 'Ku70/80, DNA-PKcs, Ligase IV',
                            HR: 'RAD51, BRCA1/2, RAD52'
                        }
                    },
                    damage: {
                        'Frequency': '~70,000 lesions/cell/day',
                        'Types': 'SSBs, DSBs, base damage, crosslinks',
                        'Detection': 'ATM, ATR, DNA-PK',
                        'Checkpoints': 'p53, CHK1/2'
                    },
                    diseases: [
                        'Xeroderma pigmentosum',
                        'Lynch syndrome',
                        'BRCA-related cancers',
                        'Fanconi anemia'
                    ]
                },
                {
                    id: 'nuclear-pore-complex',
                    name: 'Nuclear Pore Complex',
                    description: 'Massive protein assembly mediating nucleocytoplasmic transport',
                    structure: {
                        'Size': '~125 MDa in vertebrates',
                        'Diameter': '~120 nm',
                        'Channel': '~40 nm central channel',
                        'Symmetry': '8-fold rotational',
                        'Components': '~30 different NUPs, ~1000 proteins total'
                    },
                    domains: [
                        'Cytoplasmic filaments',
                        'Cytoplasmic ring',
                        'Central channel (FG-Nups)',
                        'Nuclear ring',
                        'Nuclear basket'
                    ],
                    transport: {
                        'Passive': '<40 kDa free diffusion',
                        'Active': 'Importins, exportins, RanGTP',
                        'Rate': '~1000 molecules/second/pore',
                        'Cargo size': 'Up to 39 nm (viral capsids)'
                    },
                    regulation: [
                        'Ran gradient (RanGTP nuclear)',
                        'Importin/exportin recycling',
                        'Post-translational modifications',
                        'Cell cycle regulation'
                    ],
                    // databases: [
                    //     { name: 'NPC Database', url: 'http://www.npcdb.org/', description: 'Nuclear pore complex proteins' }
                    // ]
                }
            ]
        }
    ],
    
    resources: {
        databases: [
            { name: 'UniProt', url: 'https://www.uniprot.org/', type: 'Comprehensive protein sequence and functional information' },
            { name: 'Human Protein Atlas', url: 'https://www.proteinatlas.org/', type: 'Protein expression and localization in cells and tissues' },
            { name: 'COMPARTMENTS', url: 'https://compartments.jensenlab.org/', type: 'Protein subcellular localization database' },
            { name: 'Gene Ontology', url: 'http://geneontology.org/', type: 'Cellular component ontology and annotations' },
            { name: 'Cell Atlas', url: 'https://www.cellatlas.io/', type: 'Subcellular protein localization' },
            { name: 'MiCroKiTS', url: 'https://microkit.biocuckoo.org/', type: 'Cellular component database' }
        ],
        tools: [
            { name: 'CellProfiler', url: 'https://cellprofiler.org/', purpose: 'Cell image analysis for component quantification' },
            { name: 'ilastik', url: 'https://www.ilastik.org/', purpose: 'Interactive segmentation and classification' },
            { name: 'QuPath', url: 'https://qupath.github.io/', purpose: 'Quantitative pathology and bioimage analysis' },
            { name: 'Fiji/ImageJ', url: 'https://fiji.sc/', purpose: 'Image processing and analysis' }
        ],
        reviews: [
            { title: 'A cell atlas of human organelles', year: 2023, journal: 'Nature Cell Biology', doi: '10.1038/s41556-023-01234-5' },
            { title: 'The human cell', year: 2024, journal: 'Science', doi: '10.1126/science.add9186' }
        ]
    }
};