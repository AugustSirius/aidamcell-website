export const cellularDynamics = {
    id: 'cellular-dynamics',
    name: 'Cellular Dynamics & Systems',
    description: 'Temporal processes, signaling networks, and emergent cellular behaviors that govern cell function, fate decisions, and responses to environmental stimuli',
    
    categories: [
        {
            id: 'signaling-networks',
            name: 'Signal Transduction',
            icon: 'fas fa-network-wired',
            color: 'primary',
            items: [
                {
                    id: 'receptor-signaling',
                    name: 'Receptor Signaling Cascades',
                    description: 'Information flow from extracellular signals through membrane receptors to cellular responses',
                    pathways: {
                        'RTK/MAPK': {
                            ligands: 'EGF, FGF, PDGF, VEGF, insulin',
                            cascade: 'RTK → Grb2/SOS → Ras → Raf → MEK → ERK',
                            outputs: 'Proliferation, differentiation, survival',
                            duration: 'Peak at 5-10 min, sustained 1-4 hours'
                        },
                        'PI3K/AKT': {
                            activation: 'RTKs, GPCRs, Ras',
                            cascade: 'PI3K → PIP3 → PDK1 → AKT',
                            targets: 'mTOR, GSK3β, FoxO, MDM2',
                            functions: 'Survival, metabolism, growth'
                        },
                        'JAK/STAT': {
                            ligands: 'Cytokines, interferons, growth hormones',
                            mechanism: 'Receptor dimerization → JAK → STAT',
                            genes: '4 JAKs, 7 STATs in mammals',
                            response: 'Direct nuclear translocation'
                        },
                        'GPCR': {
                            receptors: '>800 in humans',
                            G_proteins: 'Gs, Gi/o, Gq/11, G12/13',
                            messengers: 'cAMP, Ca²⁺, DAG, IP3',
                            regulation: 'β-arrestins, GRKs, RGS proteins'
                        },
                        'TGF-β/SMAD': {
                            ligands: 'TGF-β, BMPs, activins',
                            receptors: 'Type I and II serine/threonine kinases',
                            cascade: 'R-SMAD → Co-SMAD → nucleus',
                            functions: 'EMT, growth inhibition, differentiation'
                        },
                        'Wnt': {
                            canonical: 'β-catenin dependent, gene transcription',
                            noncanonical: 'PCP pathway, Ca²⁺ pathway',
                            components: 'Frizzled, LRP5/6, Dishevelled',
                            outputs: 'Cell fate, proliferation, polarity'
                        },
                        'Notch': {
                            mechanism: 'Juxtacrine, proteolytic cleavage',
                            ligands: 'Delta, Jagged/Serrate',
                            process: 'NICD release → nuclear translocation',
                            functions: 'Lateral inhibition, boundary formation'
                        }
                    },
                    dynamics: {
                        'Activation': '10 seconds - 2 minutes',
                        'Peak response': '5-30 minutes typical',
                        'Adaptation': '30 minutes - 4 hours',
                        'Transcriptional changes': '30 minutes - 8 hours',
                        'Protein expression': '2-24 hours'
                    },
                    emergentProperties: [
                        'Ultrasensitivity (Hill coefficient >1)',
                        'Bistability (toggle switches)',
                        'Oscillations (NF-κB, p53, ERK)',
                        'Perfect adaptation (robust to perturbations)',
                        'Signal amplification (10⁴-10⁶ fold)'
                    ],
                    crosstalk: {
                        'Shared components': 'Average 3-5 pathways per protein',
                        'Integration points': 'Ras, PI3K, calcium',
                        'Feedback loops': 'Negative (adaptation), positive (commitment)',
                        'Spatial organization': 'Scaffolds, compartments'
                    },
                    databases: [
                        { name: 'KEGG Pathways', url: 'https://www.kegg.jp/kegg/pathway.html', description: 'Comprehensive pathway maps' },
                        { name: 'Reactome', url: 'https://reactome.org/', description: 'Pathway database and visualization' },
                        { name: 'WikiPathways', url: 'https://www.wikipathways.org/', description: 'Community-curated pathways' },
                        { name: 'PathwayCommons', url: 'https://www.pathwaycommons.org/', description: 'Integrated pathway data' }
                    ]
                },
                {
                    id: 'calcium-signaling',
                    name: 'Calcium Dynamics',
                    description: 'Universal second messenger system controlling diverse cellular processes',
                    mechanisms: [
                        'Store-operated calcium entry (SOCE) via STIM/Orai',
                        'Voltage-gated calcium channels (L, N, P/Q, R, T types)',
                        'IP3 receptors - ER calcium release',
                        'Ryanodine receptors - CICR mechanism',
                        'TRP channels - diverse stimuli sensing',
                        'Mitochondrial calcium uniporter (MCU)'
                    ],
                    patterns: {
                        'Baseline': '50-100 nM cytosolic',
                        'Activated': '500-1000 nM peaks',
                        'Oscillations': '0.002-2 Hz frequency',
                        'Waves': '10-100 µm/s propagation',
                        'Sparks': 'Local 2-20 µm, 10-200 ms events',
                        'Plateaus': 'Sustained elevation for minutes'
                    },
                    encoding: {
                        'Amplitude modulation': 'Strength of response',
                        'Frequency modulation': 'Different gene programs',
                        'Spatial patterns': 'Local vs global signals',
                        'Duration': 'Transient vs sustained'
                    },
                    effectors: [
                        'Calmodulin - primary Ca²⁺ sensor',
                        'CaMKII - memory and plasticity',
                        'Calcineurin - NFAT activation',
                        'PKC - DAG synergy',
                        'Calpains - proteolysis'
                    ],
                    cellularFunctions: [
                        'Muscle contraction (excitation-contraction)',
                        'Neurotransmitter release (synaptic vesicles)',
                        'Gene transcription (CREB, NFAT, MEF2)',
                        'Cell migration (focal adhesions)',
                        'Fertilization (egg activation)',
                        'Apoptosis (mitochondrial pathway)'
                    ],
                    buffering: {
                        'Proteins': 'Parvalbumin, calbindin, calretinin',
                        'ER': '~500 µM storage capacity',
                        'Mitochondria': 'Rapid uptake, slow release',
                        'Extrusion': 'PMCA, NCX (~200 Ca²⁺/s)'
                    },
                    databases: [
                        { name: 'Guide to Pharmacology', url: 'https://www.guidetopharmacology.org/', description: 'Calcium channel database' }
                    ]
                },
                {
                    id: 'metabolic-signaling',
                    name: 'Metabolic Signaling',
                    description: 'Integration of nutrient availability with cellular growth and survival decisions',
                    sensors: {
                        'AMPK': {
                            activation: 'AMP/ATP ratio >10',
                            targets: 'ACC, mTOR, PGC-1α, p53',
                            effects: 'ATP production ↑, consumption ↓',
                            drugs: 'Metformin, AICAR'
                        },
                        'mTOR': {
                            complexes: 'mTORC1 (growth), mTORC2 (survival)',
                            inputs: 'Amino acids, growth factors, energy',
                            outputs: 'Translation, lipogenesis, autophagy ⊣',
                            inhibitors: 'Rapamycin, everolimus'
                        },
                        'HIF': {
                            regulation: 'O₂-dependent hydroxylation',
                            targets: 'VEGF, GLUT1, glycolytic enzymes',
                            threshold: '<5% O₂ activation',
                            diseases: 'Cancer, ischemia'
                        },
                        'Sirtuins': {
                            cofactor: 'NAD⁺ dependent',
                            substrates: 'Histones, metabolic enzymes',
                            functions: 'Deacetylation, longevity',
                            number: 'SIRT1-7 in mammals'
                        },
                        'PKA': {
                            activation: 'cAMP from adenylyl cyclase',
                            fasting: 'Glucagon → glycogenolysis',
                            fed: 'Insulin opposes PKA',
                            compartments: 'AKAPs localize signaling'
                        }
                    },
                    metabolicStates: [
                        'Glycolytic (Warburg) - cancer, activation',
                        'OXPHOS - differentiated cells, quiescence',
                        'FAO - memory T cells, M2 macrophages',
                        'Glutaminolysis - proliferating cells',
                        'One-carbon - nucleotide synthesis'
                    ],
                    integration: {
                        'Fed state': 'Insulin → anabolism',
                        'Fasted state': 'Glucagon → catabolism',
                        'Exercise': 'AMPK → substrate mobilization',
                        'Circadian': '~50% metabolic genes oscillate'
                    },
                    diseases: [
                        'Diabetes - insulin resistance',
                        'Cancer - metabolic reprogramming',
                        'Neurodegeneration - mitochondrial dysfunction',
                        'NAFLD - lipid accumulation'
                    ],
                    databases: [
                        { name: 'HMDB', url: 'https://hmdb.ca/', description: 'Human Metabolome Database' },
                        { name: 'MetaboAnalyst', url: 'https://www.metaboanalyst.ca/', description: 'Metabolomics analysis' }
                    ]
                },
                {
                    id: 'stress-signaling',
                    name: 'Stress Response Pathways',
                    description: 'Cellular mechanisms for detecting and responding to various stressors',
                    types: {
                        'ER stress': {
                            sensors: 'IRE1, PERK, ATF6',
                            response: 'UPR - protein folding capacity ↑',
                            outcomes: 'Adaptation, autophagy, apoptosis',
                            markers: 'BiP, CHOP, XBP1s'
                        },
                        'Oxidative stress': {
                            sensors: 'KEAP1/NRF2 system',
                            ROS_sources: 'Mitochondria, NOX, peroxisomes',
                            antioxidants: 'GSH, SOD, catalase, peroxiredoxins',
                            damage: 'Lipid peroxidation, DNA breaks, protein carbonylation'
                        },
                        'DNA damage': {
                            sensors: 'ATM (DSBs), ATR (SSBs, stalled forks)',
                            checkpoints: 'G1/S, intra-S, G2/M',
                            effectors: 'p53, CHK1/2, H2AX',
                            repair: 'HR, NHEJ, BER, NER, MMR'
                        },
                        'Heat shock': {
                            sensor: 'HSF1 trimerization',
                            response: 'HSP expression ↑',
                            proteins: 'Hsp70, Hsp90, small HSPs',
                            protection: 'Protein refolding, aggregation prevention'
                        },
                        'Hyperosmotic': {
                            sensors: 'WNK kinases',
                            response: 'p38/JNK activation',
                            adaptation: 'Organic osmolyte accumulation',
                            genes: 'AR, BGT1, SMIT'
                        }
                    },
                    integration: 'ISR (Integrated Stress Response) via eIF2α phosphorylation',
                    databases: [
                        { name: 'Stressrome', url: 'http://stressrome.org/', description: 'Stress response database' }
                    ]
                }
            ]
        },
        {
            id: 'cell-cycle-dynamics',
            name: 'Cell Cycle & Division',
            icon: 'fas fa-circle-notch',
            color: 'secondary',
            items: [
                {
                    id: 'cell-cycle-control',
                    name: 'Cell Cycle Regulation',
                    description: 'Ordered series of events leading to cell division with multiple checkpoints ensuring genomic integrity',
                    phases: {
                        'G1': {
                            duration: '10-12 hours (variable)',
                            events: 'Growth, organelle duplication',
                            cyclins: 'Cyclin D + CDK4/6',
                            checkpoint: 'Restriction point (R)',
                            commitment: 'Rb phosphorylation'
                        },
                        'S': {
                            duration: '6-8 hours',
                            events: 'DNA replication',
                            cyclins: 'Cyclin E/A + CDK2',
                            regulation: 'Origin licensing and firing',
                            quality: 'Replication checkpoint'
                        },
                        'G2': {
                            duration: '3-4 hours',
                            events: 'Preparation for mitosis',
                            cyclins: 'Cyclin A/B + CDK1',
                            checkpoint: 'DNA damage check',
                            centrosomes: 'Duplication completion'
                        },
                        'M': {
                            duration: '~1 hour',
                            events: 'Nuclear division, cytokinesis',
                            cyclins: 'Cyclin B + CDK1',
                            checkpoint: 'Spindle assembly (SAC)',
                            APC_C: 'Anaphase promotion'
                        },
                        'G0': {
                            state: 'Quiescent/differentiated',
                            cells: 'Most adult cells',
                            reentry: 'Growth signals required',
                            markers: 'p27 high, Ki67 negative'
                        }
                    },
                    checkpoints: [
                        'G1/S - DNA damage, cell size, nutrients',
                        'Intra-S - Replication fork stalling',
                        'G2/M - DNA damage, unreplicated DNA',
                        'Spindle - All kinetochores attached'
                    ],
                    molecularOscillators: {
                        'CDK activity': 'Low(G1) → Med(S) → High(G2/M)',
                        'Cyclins': 'D(G1) → E(G1/S) → A(S/G2) → B(M)',
                        'CKIs': 'p21(Cip1), p27(Kip1), p57(Kip2)',
                        'Pocket proteins': 'Rb, p107, p130',
                        'APC/C': 'Cdh1(G1), Cdc20(M)',
                        'Phosphatases': 'CDC25A/B/C, PP2A'
                    },
                    singleCellVariability: {
                        'G1 duration': 'CV ~30-40%',
                        'Immediate vs delayed': 'Bifurcation in daughter cells',
                        'p21 dynamics': 'Pulsatile after DNA damage',
                        'Mitogen memory': 'Integration over time'
                    },
                    diseases: [
                        'Cancer - checkpoint bypass',
                        'Developmental disorders',
                        'Neurodegenerative - aberrant re-entry'
                    ],
                    databases: [
                        { name: 'Cyclebase', url: 'https://cyclebase.org/', description: 'Cell cycle gene database' },
                        { name: 'mitocheck', url: 'https://www.mitocheck.org/', description: 'Cell division genes' }
                    ]
                },
                {
                    id: 'mitosis-mechanics',
                    name: 'Mitotic Machinery',
                    description: 'Mechanical and molecular processes ensuring accurate chromosome segregation',
                    stages: [
                        'Prophase - Chromatin condensation (condensin), centrosome migration',
                        'Prometaphase - NEB, kinetochore capture, chromosome congression',
                        'Metaphase - Bi-orientation, alignment at equator, SAC satisfaction',
                        'Anaphase A - Sister separation (separase), poleward movement',
                        'Anaphase B - Spindle elongation, pole separation',
                        'Telophase - Nuclear reformation, chromatin decondensation',
                        'Cytokinesis - Contractile ring, abscission'
                    ],
                    spindle: {
                        'Components': 'Tubulin, motor proteins, MAPs',
                        'Types': 'k-fibers, interpolar, astral MTs',
                        'Assembly': 'Search-and-capture, RanGTP gradient',
                        'Forces': 'Pushing (kinesin-5), pulling (dynein)',
                        'Length': '10-25 µm typical'
                    },
                    kinetochore: {
                        'Structure': '>100 proteins, ~1 MDa',
                        'Layers': 'Inner (CENP-A) → outer (Ndc80)',
                        'Attachment': '20-30 MTs per kinetochore',
                        'Tension': '~100 pN per kinetochore',
                        'Error correction': 'Aurora B gradient'
                    },
                    forces: {
                        'Spindle tension': '100-700 pN total',
                        'Cortical tension': '1-5 nN/µm',
                        'Chromosome velocity': '0.5-2 µm/min',
                        'Cytokinetic ring': '10-100 nN total force'
                    },
                    regulation: {
                        'SAC proteins': 'Mad2, BubR1, Bub1, Mps1',
                        'APC/C-Cdc20': 'Ubiquitinates cyclin B, securin',
                        'Separase': 'Cleaves cohesin',
                        'Aurora/Polo kinases': 'Multiple mitotic roles'
                    },
                    errors: {
                        'Rate': '1/100 to 1/1000 divisions',
                        'Types': 'Lagging, bridges, merotelic',
                        'Consequences': 'Aneuploidy, micronuclei',
                        'CIN': 'Chromosomal instability in cancer'
                    },
                    asymmetry: {
                        'Mechanisms': 'Spindle positioning, fate determinants',
                        'Examples': 'Stem cells, neuroblasts, budding yeast',
                        'Regulators': 'aPKC, Par proteins, Inscuteable'
                    },
                    databases: [
                        { name: 'MitoCheck', url: 'https://www.mitocheck.org/', description: 'Mitosis gene functions' }
                    ]
                },
                {
                    id: 'dna-replication',
                    name: 'DNA Replication Dynamics',
                    description: 'Semi-conservative duplication of the genome with high fidelity',
                    machinery: {
                        'Helicase': 'CMG complex, 10-50 bp/s unwinding',
                        'Primase': 'Pol α-primase, RNA primers',
                        'Polymerases': 'Pol ε (leading), Pol δ (lagging)',
                        'Clamp': 'PCNA processivity factor',
                        'Ligase': 'DNA ligase I, Okazaki joining'
                    },
                    statistics: {
                        'Speed': '50 nucleotides/second/fork',
                        'Origins': '30,000-50,000 licensed, 10-20% fired',
                        'Fork progression': '0.5-2 kb/min',
                        'Error rate': '10⁻⁹-10⁻¹⁰ per base',
                        'Processivity': '>100 kb without dissociation'
                    },
                    regulation: [
                        'Origin licensing (G1) - MCM2-7 loading',
                        'Origin firing (S) - DDK, CDK activation',
                        'Replication timing - Early (active), late (heterochromatin)',
                        'Fork protection - RPA, RAD51',
                        'Checkpoint activation - ATR/CHK1'
                    ],
                    challenges: {
                        'Conflicts': 'Transcription-replication collisions',
                        'Obstacles': 'DNA damage, secondary structures',
                        'Difficult regions': 'Repetitive DNA, telomeres',
                        'Fork stalling': 'dNTP depletion, DNA damage'
                    },
                    stress: {
                        'Causes': 'Oncogenes, nucleotide imbalance',
                        'Response': 'Fork reversal, dormant origin firing',
                        'Resolution': 'HR repair, break-induced replication',
                        'Consequences': 'Genomic instability → cancer'
                    },
                    databases: [
                        { name: 'Replication Domain', url: 'https://www.replicationdomain.com/', description: 'Replication timing data' }
                    ]
                },
                {
                    id: 'cytokinesis',
                    name: 'Cytokinesis',
                    description: 'Physical division of the cytoplasm to produce two daughter cells',
                    mechanisms: {
                        'Animals': 'Contractile ring (actomyosin)',
                        'Plants': 'Cell plate (phragmoplast)',
                        'Fungi': 'Septum formation',
                        'Bacteria': 'Z-ring (FtsZ)'
                    },
                    contractileRing: {
                        'Components': 'Actin, myosin II, septins',
                        'Assembly': 'Anillin, formins, RhoA',
                        'Contraction': '~0.5 µm/min rate',
                        'Force': '10-100 nN total'
                    },
                    positioning: {
                        'Signals': 'Spindle midzone, astral MTs',
                        'Proteins': 'CPC (Aurora B), centralspindlin',
                        'Symmetry': 'Usually equal division',
                        'Asymmetry': 'Polar bodies, stem cells'
                    },
                    abscission: {
                        'Bridge': 'Midbody, 1-2 µm diameter',
                        'Duration': '1-3 hours after furrow',
                        'ESCRT': 'Membrane scission complex',
                        'Checkpoint': 'Aurora B delays if lagging DNA'
                    }
                }
            ]
        },
        {
            id: 'cellular-mechanics',
            name: 'Cellular Mechanics',
            icon: 'fas fa-cogs',
            color: 'success',
            items: [
                {
                    id: 'cytoskeletal-dynamics',
                    name: 'Cytoskeletal Dynamics',
                    description: 'Dynamic polymer networks generating forces and maintaining cell structure',
                    systems: {
                        'Actin': {
                            polymerization: {
                                'Rate': '11.6 µM⁻¹s⁻¹ (barbed end)',
                                'Critical concentration': '0.1 µM (barbed), 0.7 µM (pointed)',
                                'ATP hydrolysis': 'Lag time ~2 seconds',
                                'Treadmilling': '0.04-0.3 µm/min'
                            },
                            structures: {
                                'Lamellipodia': 'Branched network, Arp2/3',
                                'Filopodia': 'Parallel bundles, fascin',
                                'Stress fibers': 'Contractile, α-actinin',
                                'Cortex': 'Spectrin, ERM proteins'
                            },
                            regulation: 'Profilin, cofilin, capping proteins, formins'
                        },
                        'Microtubules': {
                            dynamics: {
                                'Growth': '15-20 µm/min',
                                'Shrinkage': '30-60 µm/min',
                                'Catastrophe': '0.005-0.05 s⁻¹',
                                'Rescue': '0.01-0.1 s⁻¹'
                            },
                            structures: {
                                'Centrosome': 'MTOC, γ-tubulin',
                                'Spindle': 'Mitotic apparatus',
                                'Cilia': 'Primary, motile',
                                'Axon': 'Parallel arrays'
                            },
                            modifications: 'Acetylation, tyrosination, polyglutamylation'
                        },
                        'Intermediate filaments': {
                            types: 'Keratins (epithelial), vimentin (mesenchymal), GFAP (glial)',
                            assembly: 'Unit-length-filament model',
                            mechanics: 'Strain-stiffening, viscoelastic',
                            turnover: 'Hours to days (much slower)'
                        }
                    },
                    forces: {
                        'Polymerization': '5-7 pN per filament',
                        'Motor proteins': '1-7 pN (single molecule)',
                        'Cell traction': '10-100 nN total',
                        'Cortical tension': '0.1-5 mN/m'
                    },
                    regulation: {
                        'Small GTPases': {
                            'RhoA': 'Stress fibers, contractility',
                            'Rac1': 'Lamellipodia, membrane ruffles',
                            'Cdc42': 'Filopodia, polarity'
                        },
                        'Kinases': 'ROCK, PAK, LIMK, MLCK',
                        'Scaffolds': 'IQGAP, WASP, WAVE'
                    },
                    mechanosensing: [
                        'Catch bonds in focal adhesions',
                        'Filamin unfolding',
                        'Spectrin unfolding in RBCs',
                        'Nuclear mechanotransduction'
                    ],
                    databases: [
                        { name: 'CytoSkeletal Signaling', url: 'https://www.cellsignal.com/pathways/cytoskeletal-signaling', description: 'Cytoskeleton pathways' }
                    ]
                },
                {
                    id: 'cell-migration',
                    name: 'Cell Migration',
                    description: 'Coordinated processes enabling directed cell movement',
                    types: {
                        'Mesenchymal': {
                            speed: '0.1-1 µm/min',
                            features: 'Strong adhesions, proteolysis',
                            cells: 'Fibroblasts, cancer invasion'
                        },
                        'Amoeboid': {
                            speed: '2-25 µm/min',
                            features: 'Weak adhesions, squeezing',
                            cells: 'Leukocytes, some cancer cells'
                        },
                        'Collective': {
                            speed: '0.01-0.5 µm/min',
                            features: 'Cell-cell junctions maintained',
                            examples: 'Epithelial sheets, neural crest'
                        },
                        'Swimming': {
                            mechanism: 'Flagella, shape changes',
                            cells: 'Sperm, some bacteria'
                        }
                    },
                    steps: [
                        'Polarization - Front-rear axis establishment',
                        'Protrusion - Actin polymerization at leading edge',
                        'Adhesion - Integrin clustering, nascent adhesions',
                        'Traction - Actomyosin contraction',
                        'De-adhesion - Rear release, FAK/Src signaling'
                    ],
                    adhesions: {
                        'Nascent': '<1 µm, <1 min lifetime',
                        'Focal complexes': '1-2 µm, few minutes',
                        'Focal adhesions': '2-5 µm, 10-60 min',
                        'Fibrillar': '>5 µm, stable',
                        'Force transmission': '5-30 kPa typical'
                    },
                    guidance: {
                        'Chemotaxis': {
                            gradients: '2-10% across cell',
                            receptors: 'GPCRs, RTKs',
                            examples: 'fMLP, chemokines, growth factors'
                        },
                        'Haptotaxis': {
                            mechanism: 'ECM gradients',
                            molecules: 'Fibronectin, laminin, collagen'
                        },
                        'Durotaxis': {
                            range: '1-100 kPa stiffness',
                            sensing: 'Focal adhesion maturation'
                        },
                        'Galvanotaxis': {
                            fields: '0.1-10 V/cm physiological',
                            mechanism: 'Voltage-gated channels, electrophoresis'
                        }
                    },
                    physiologicalRoles: [
                        'Embryonic development - gastrulation, organogenesis',
                        'Immune surveillance - leukocyte trafficking',
                        'Wound healing - re-epithelialization',
                        'Cancer metastasis - invasion, intravasation'
                    ],
                    databases: [
                        { name: 'Cell Migration Gateway', url: 'https://www.cellmigration.org/', description: 'Migration resources' }
                    ]
                },
                {
                    id: 'mechanotransduction',
                    name: 'Mechanotransduction',
                    description: 'Conversion of mechanical stimuli into biochemical signals',
                    sensors: [
                        'Focal adhesions - Integrin clustering, talin unfolding',
                        'Stretch-activated channels - Piezo1/2, TRP channels',
                        'Cytoskeleton - Spectrin unfolding, catch bonds',
                        'Nuclear envelope - Lamin A/C, LINC complex',
                        'Primary cilium - Bending → Ca²⁺ influx',
                        'Glycocalyx - Shear stress sensing'
                    ],
                    mechanisms: {
                        'Conformational changes': 'Protein unfolding reveals sites',
                        'Clustering': 'Receptor aggregation',
                        'Membrane tension': 'Lipid packing, channel gating',
                        'Cytoskeletal prestress': 'Tensegrity model'
                    },
                    responses: {
                        'Immediate (<1s)': 'Ion channel opening, Ca²⁺ flux',
                        'Fast (seconds)': 'Protein phosphorylation, unfolding',
                        'Intermediate (minutes)': 'Cytoskeletal remodeling',
                        'Slow (hours)': 'Gene expression, differentiation'
                    },
                    forceSensitivity: {
                        'Single molecules': '1-10 pN unfolding',
                        'Focal adhesions': '5-30 kPa substrate',
                        'Cell deformation': '1-10% strain',
                        'Shear stress': '1-30 dyn/cm² (vessels)'
                    },
                    stiffnessRange: {
                        'Brain': '0.1-1 kPa',
                        'Fat': '2-4 kPa',
                        'Muscle': '10-20 kPa',
                        'Bone': '25-40 GPa'
                    },
                    pathways: [
                        'YAP/TAZ - Hippo pathway mechanosensing',
                        'MRTF-A - Actin-regulated transcription',
                        'NF-κB - Shear stress response',
                        'β-catenin - E-cadherin tension sensing'
                    ],
                    diseases: [
                        'Atherosclerosis - disturbed flow',
                        'Fibrosis - excessive ECM stiffening',
                        'Cancer - altered mechanosensing',
                        'Muscular dystrophies - membrane fragility'
                    ],
                    databases: [
                        { name: 'MechanoBase', url: 'https://mechanobase.org/', description: 'Mechanobiology database' }
                    ]
                },
                {
                    id: 'cell-shape',
                    name: 'Cell Shape & Polarity',
                    description: 'Mechanisms controlling cell morphology and asymmetry',
                    determinants: [
                        'Cortical tension - actomyosin contractility',
                        'Adhesion geometry - ECM patterning',
                        'Membrane trafficking - exo/endocytosis balance',
                        'Osmotic pressure - water/ion regulation'
                    ],
                    polarity: {
                        'Apical-basal': {
                            complexes: 'Par, Crumbs, Scribble',
                            cells: 'Epithelia, neurons'
                        },
                        'Planar': {
                            pathway: 'PCP (Frizzled, Vangl)',
                            examples: 'Hair cells, wing cells'
                        },
                        'Front-rear': {
                            regulators: 'PI3K-PTEN, Cdc42',
                            context: 'Migration, chemotaxis'
                        }
                    },
                    shapeMaintenance: {
                        'Red blood cells': 'Spectrin network, biconcave',
                        'Neurons': 'Axon initial segment, dendrites',
                        'Epithelia': 'Cell-cell junctions, columnar/squamous'
                    }
                }
            ]
        },
        {
            id: 'cellular-communication',
            name: 'Cell-Cell Communication',
            icon: 'fas fa-comments',
            color: 'warning',
            items: [
                {
                    id: 'direct-contact',
                    name: 'Direct Cell Contact',
                    description: 'Communication through physical interactions between adjacent cells',
                    mechanisms: [
                        'Gap junctions - Direct cytoplasmic coupling',
                        'Tunneling nanotubes - Long-range connections',
                        'Cell adhesion molecules - Signal transduction',
                        'Notch-Delta - Juxtacrine signaling',
                        'Ephrin-Eph - Bidirectional signaling',
                        'Immune synapses - T cell activation'
                    ],
                    gapJunctions: {
                        'Structure': 'Connexin hexamers (connexons)',
                        'Pore size': '1.5 nm diameter',
                        'Permeability': '<1 kDa molecules',
                        'Conductance': '50-300 pS',
                        'Regulation': 'pH, Ca²⁺, voltage, phosphorylation',
                        'Functions': 'Metabolic coupling, electrical synapses'
                    },
                    tunnelingNanotubes: {
                        'Diameter': '50-200 nm',
                        'Length': 'Up to 100 µm',
                        'Contents': 'Organelles, proteins, RNA',
                        'Formation': 'Actin-driven protrusion',
                        'Functions': 'Rescue, pathogen spread'
                    },
                    notchSignaling: {
                        'Ligands': 'Delta1-4, Jagged1-2',
                        'Receptors': 'Notch1-4',
                        'Mechanism': 'Force-dependent unfolding',
                        'Output': 'NICD → transcription',
                        'Functions': 'Lateral inhibition, boundaries'
                    },
                    immuneSynapse: {
                        'Structure': 'cSMAC, pSMAC, dSMAC',
                        'Duration': '30 min - hours',
                        'Molecules': 'TCR, MHC, costimulatory',
                        'Outcome': 'T cell activation/tolerance'
                    },
                    databases: [
                        { name: 'CellPhoneDB', url: 'https://www.cellphonedb.org/', description: 'Cell-cell interaction database' }
                    ]
                },
                {
                    id: 'paracrine-signaling',
                    name: 'Paracrine Signaling',
                    description: 'Local signaling between nearby cells through secreted factors',
                    molecules: {
                        'Growth factors': {
                            families: 'FGF (23), VEGF (5), PDGF (4), EGF',
                            range: '20-200 µm typical',
                            regulation: 'ECM binding, proteolysis'
                        },
                        'Cytokines': {
                            types: 'Interleukins (40), interferons, TNF',
                            functions: 'Inflammation, immunity',
                            pleiotropy: 'Multiple cell targets'
                        },
                        'Chemokines': {
                            families: 'CC, CXC, CX3C, XC',
                            gradients: 'Haptotaxis, chemotaxis',
                            receptors: 'GPCRs (20 types)'
                        },
                        'Morphogens': {
                            examples: 'BMP, Shh, Wnt, retinoic acid',
                            gradients: 'French flag model',
                            range: '100-500 µm'
                        }
                    },
                    diffusion: {
                        'Coefficient': '10-100 µm²/s in tissue',
                        'Range': 'λ = √(D/k) decay length',
                        'Barriers': 'ECM, binding proteins',
                        'Transport': 'Transcytosis, cytonemes'
                    },
                    regulation: [
                        'ECM sequestration - HSPGs, fibronectin',
                        'Proteolytic activation - MMPs, ADAMs',
                        'Decoy receptors - Soluble forms',
                        'Feedback loops - Auto/paracrine'
                    ],
                    ECMinteractions: {
                        'Binding': 'Heparan sulfate, fibronectin',
                        'Presentation': 'Concentrated, oriented',
                        'Release': 'Protease cleavage',
                        'Modulation': 'Bioavailability control'
                    },
                    databases: [
                        { name: 'MatrisomeDB', url: 'http://matrisomeproject.mit.edu/', description: 'ECM protein database' }
                    ]
                },
                {
                    id: 'extracellular-vesicles',
                    name: 'Extracellular Vesicles',
                    description: 'Membrane-bound particles mediating intercellular communication',
                    types: {
                        'Exosomes': {
                            size: '30-150 nm',
                            origin: 'Endosomal (MVB fusion)',
                            markers: 'CD9, CD63, CD81, TSG101',
                            release: 'Constitutive + regulated'
                        },
                        'Microvesicles': {
                            size: '100-1000 nm',
                            origin: 'Plasma membrane budding',
                            triggers: 'Ca²⁺, cell activation',
                            enrichment: 'PS externalization'
                        },
                        'Apoptotic bodies': {
                            size: '1-5 µm',
                            origin: 'Dying cells',
                            content: 'Nuclear fragments',
                            clearance: 'Phagocytosis'
                        },
                        'Large oncosomes': {
                            size: '1-10 µm',
                            cells: 'Cancer cells',
                            mechanism: 'Membrane blebbing'
                        }
                    },
                    cargo: [
                        'Proteins - 1000s of types',
                        'mRNA - functional transfer',
                        'miRNA - gene regulation',
                        'Long ncRNA - regulatory',
                        'DNA - genomic, mitochondrial',
                        'Lipids - bioactive mediators'
                    ],
                    biogenesis: {
                        'ESCRT pathway': 'MVB formation',
                        'ESCRT-independent': 'Ceramide, tetraspanins',
                        'Loading': 'KFERQ motifs, SUMOylation',
                        'Release triggers': 'Ca²⁺, pH, hypoxia'
                    },
                    uptake: {
                        'Mechanisms': 'Endocytosis, fusion, phagocytosis',
                        'Specificity': 'Receptor-ligand, tissue tropism',
                        'Barriers': 'Glycocalyx, immune clearance'
                    },
                    functions: [
                        'Waste removal - protein aggregates',
                        'Immune regulation - antigen presentation',
                        'Development - morphogen transport',
                        'Cancer - metastatic niche preparation',
                        'Regeneration - stem cell activation'
                    ],
                    concentration: '10¹⁰-10¹² vesicles/mL in blood',
                    databases: [
                        { name: 'ExoCarta', url: 'http://www.exocarta.org/', description: 'Exosome database' },
                        { name: 'Vesiclepedia', url: 'http://www.microvesicles.org/', description: 'Extracellular vesicle database' },
                        { name: 'EVpedia', url: 'http://evpedia.info/', description: 'EV proteins, RNAs, lipids' }
                    ]
                },
                {
                    id: 'synaptic-transmission',
                    name: 'Synaptic Transmission',
                    description: 'Specialized cell-cell communication in the nervous system',
                    types: {
                        'Chemical': {
                            neurotransmitters: 'Glutamate, GABA, dopamine, etc.',
                            delay: '0.5-5 ms',
                            plasticity: 'LTP, LTD',
                            modulation: 'Presynaptic, postsynaptic'
                        },
                        'Electrical': {
                            structure: 'Gap junctions',
                            speed: 'No synaptic delay',
                            bidirectional: 'Usually',
                            locations: 'Retina, heart, glia'
                        }
                    },
                    vesicleRelease: {
                        'Pools': 'RRP (~10), recycling (~100), reserve',
                        'Fusion': 'SNAREs, synaptotagmin',
                        'Modes': 'Synchronous, asynchronous, spontaneous',
                        'Probability': '0.1-0.9 per action potential'
                    },
                    plasticity: {
                        'Short-term': 'Facilitation, depression (ms-s)',
                        'Long-term': 'LTP, LTD (hours-lifetime)',
                        'Structural': 'Spine formation/elimination',
                        'Homeostatic': 'Scaling, metaplasticity'
                    }
                }
            ]
        },
        {
            id: 'cellular-states',
            name: 'Cellular States & Transitions',
            icon: 'fas fa-exchange-alt',
            color: 'info',
            items: [
                {
                    id: 'differentiation',
                    name: 'Cell Differentiation',
                    description: 'Progressive restriction of developmental potential and acquisition of specialized functions',
                    mechanisms: [
                        'Master TF activation - MyoD, Sox2, Oct4',
                        'Chromatin remodeling - BAF complex switching',
                        'DNA methylation - CpG islands, enhancers',
                        'Histone modifications - H3K27me3 → H3K4me3',
                        'Signal integration - Morphogen gradients',
                        'Feedback loops - Auto-activation, mutual inhibition'
                    ],
                    waddingtonLandscape: {
                        'Stem cells': 'High potential energy peak',
                        'Progenitors': 'Bifurcation points, valleys',
                        'Differentiated': 'Deep stable valleys',
                        'Transdifferentiation': 'Valley-to-valley crossing',
                        'Dedifferentiation': 'Climbing back up'
                    },
                    hierarchies: {
                        'Hematopoiesis': {
                            stem: 'HSC → MPP',
                            progenitors: 'CMP, CLP, MEP',
                            mature: '10+ blood cell types'
                        },
                        'Neurogenesis': {
                            stem: 'NSC → NPC',
                            fate: 'Neurons, astrocytes, oligodendrocytes',
                            diversity: '100s of neuron subtypes'
                        },
                        'Myogenesis': {
                            master: 'MyoD, Myf5',
                            progression: 'Myoblast → myotube → fiber'
                        }
                    },
                    reversibility: {
                        'iPSCs': 'Oct4, Sox2, Klf4, c-Myc (OSKM)',
                        'Efficiency': '0.01-1% typically',
                        'Partial': 'Rejuvenation without dedifferentiation',
                        'Direct': 'Transdifferentiation between lineages'
                    },
                    singleCellInsights: [
                        'Continuous trajectories vs discrete states',
                        'Multilineage priming',
                        'Stochastic vs deterministic',
                        'Cell-cell variability in timing'
                    ],
                    databases: [
                        { name: 'CellNet', url: 'http://cellnet.hms.harvard.edu/', description: 'Cell fate engineering' },
                        { name: 'ESCAPE', url: 'https://www.maayanlab.net/ESCAPE/', description: 'Stem cell signatures' }
                    ]
                },
                {
                    id: 'cell-death',
                    name: 'Programmed Cell Death',
                    description: 'Regulated cellular demise with distinct morphological and biochemical features',
                    types: {
                        'Apoptosis': {
                            pathways: 'Intrinsic (mitochondrial), extrinsic (death receptor)',
                            morphology: 'Shrinkage, blebbing, fragmentation',
                            biochemistry: 'Caspases, PS exposure, DNA laddering',
                            clearance: 'Efferocytosis, "find me/eat me" signals',
                            duration: '2-4 hours decision, 30-60 min execution'
                        },
                        'Necroptosis': {
                            triggers: 'TNF + caspase inhibition',
                            mediators: 'RIPK1, RIPK3, MLKL',
                            features: 'Membrane rupture, DAMPs release',
                            inflammation: 'Strong pro-inflammatory'
                        },
                        'Pyroptosis': {
                            triggers: 'Inflammasome activation',
                            executioner: 'Gasdermin pores',
                            release: 'IL-1β, IL-18',
                            cells: 'Immune cells primarily'
                        },
                        'Ferroptosis': {
                            mechanism: 'Lipid peroxidation',
                            requirements: 'Iron, PUFAs',
                            inhibition: 'GPX4, vitamin E',
                            morphology: 'Mitochondrial shrinkage'
                        },
                        'Autophagy-dependent': {
                            process: 'Excessive self-digestion',
                            regulation: 'mTOR, Beclin-1',
                            context: 'Starvation, stress'
                        },
                        'Entosis': {
                            mechanism: 'Cell-in-cell death',
                            context: 'Matrix detachment',
                            outcome: 'Inner cell death'
                        }
                    },
                    apoptosisRegulation: {
                        'BCL-2 family': {
                            antiapoptotic: 'BCL-2, BCL-XL, MCL-1',
                            proapoptotic: 'BAX, BAK',
                            BH3only: 'BID, BIM, PUMA, NOXA'
                        },
                        'Caspases': {
                            initiator: 'Caspase-8, -9, -10',
                            executioner: 'Caspase-3, -6, -7',
                            inflammatory: 'Caspase-1, -4, -5, -11'
                        },
                        'IAPs': 'XIAP, cIAP1/2, survivin',
                        'Death receptors': 'Fas, TNFR1, DR4/5'
                    },
                    detection: {
                        'Early': 'PS exposure (Annexin V)',
                        'Nuclear': 'TUNEL, γH2AX',
                        'Membrane': 'PI uptake',
                        'Mitochondrial': 'ΔΨm loss, cytochrome c'
                    },
                    physiologicalRoles: [
                        'Development - digit separation, neuron pruning',
                        'Homeostasis - cell turnover',
                        'Immunity - negative selection, infection control',
                        'Cancer suppression - damaged cell removal'
                    ],
                    databases: [
                        { name: 'Deathbase', url: 'http://www.deathbase.org/', description: 'Cell death database' }
                    ]
                },
                {
                    id: 'stress-responses',
                    name: 'Cellular Stress Responses',
                    description: 'Adaptive programs activated by various cellular stresses',
                    stressTypes: {
                        'Heat shock': {
                            sensor: 'HSF1 trimerization',
                            response: 'HSP expression',
                            threshold: '>42°C typically',
                            recovery: 'Thermotolerance'
                        },
                        'ER stress': {
                            sensors: 'IRE1, PERK, ATF6',
                            response: 'UPR activation',
                            outcomes: 'Folding↑, translation↓, ERAD↑',
                            chronic: 'CHOP → apoptosis'
                        },
                        'Oxidative': {
                            sensor: 'KEAP1/NRF2',
                            response: 'Antioxidant genes',
                            targets: 'GSH, NQO1, HO-1',
                            damage: '8-oxoG, protein carbonyls'
                        },
                        'DNA damage': {
                            sensors: 'ATM, ATR, DNA-PK',
                            checkpoint: 'p53 activation',
                            outcomes: 'Repair, arrest, death',
                            threshold: '~20 DSBs → death'
                        },
                        'Hypoxia': {
                            sensor: 'PHD → HIF stabilization',
                            threshold: '<5% O₂',
                            targets: 'Glycolysis, angiogenesis',
                            adaptation: 'Metabolic reprogramming'
                        },
                        'Nutrient': {
                            glucose: 'AMPK activation',
                            amino_acids: 'GCN2, mTOR inhibition',
                            response: 'Autophagy, metabolism shift'
                        }
                    },
                    integratedStressResponse: {
                        'Kinases': 'PERK, PKR, GCN2, HRI',
                        'Target': 'eIF2α phosphorylation',
                        'Effect': 'Global translation ↓',
                        'Selective': 'ATF4, CHOP translation ↑'
                    },
                    outcomes: [
                        'Adaptation - Stress protein expression',
                        'Autophagy - Cellular cleaning',
                        'Senescence - Permanent arrest',
                        'Death - If overwhelming',
                        'Transformation - Oncogenic stress'
                    ],
                    hormesis: {
                        'Concept': 'Low stress → enhanced resilience',
                        'Examples': 'Exercise, caloric restriction',
                        'Mechanisms': 'Mitohormesis, xenohormesis'
                    },
                    databases: [
                        { name: 'StressDB', url: 'https://stressdb.org/', description: 'Stress response database' }
                    ]
                },
                {
                    id: 'senescence',
                    name: 'Cellular Senescence',
                    description: 'Stable cell cycle arrest with altered secretory phenotype',
                    triggers: [
                        'Replicative - Telomere attrition',
                        'Oncogene-induced - Ras, BRAF',
                        'DNA damage - Persistent DDR',
                        'Oxidative stress - ROS accumulation',
                        'Mitochondrial dysfunction'
                    ],
                    markers: {
                        'Growth arrest': 'p16, p21, p53',
                        'SA-β-gal': 'pH 6.0 activity',
                        'Morphology': 'Enlarged, flattened',
                        'Heterochromatin': 'SAHF formation',
                        'DNA damage': 'γH2AX foci'
                    },
                    SASP: {
                        'Components': 'IL-6, IL-8, MMPs, growth factors',
                        'Regulation': 'NF-κB, C/EBPβ, mTOR',
                        'Functions': 'Paracrine senescence, immune recruitment',
                        'Timeline': 'Develops over days-weeks'
                    },
                    roles: {
                        'Beneficial': 'Tumor suppression, wound healing',
                        'Detrimental': 'Aging, inflammation, fibrosis',
                        'Development': 'Programmed senescence'
                    }
                }
            ]
        }
    ],
    
    resources: {
        databases: [
            { name: 'Reactome', url: 'https://reactome.org/', type: 'Comprehensive pathway database' },
            { name: 'KEGG', url: 'https://www.kegg.jp/', type: 'Pathways and systems biology' },
            { name: 'WikiPathways', url: 'https://www.wikipathways.org/', type: 'Community-curated pathways' },
            { name: 'SignaLink', url: 'http://signalink.org/', type: 'Signaling network resource' },
            { name: 'BioGRID', url: 'https://thebiogrid.org/', type: 'Protein and genetic interactions' },
            { name: 'STRING', url: 'https://string-db.org/', type: 'Protein-protein interactions' },
            { name: 'PhosphoSitePlus', url: 'https://www.phosphosite.org/', type: 'PTMs and signaling' }
        ],
        tools: [
            { name: 'CellNOpt', url: 'https://saezlab.github.io/CellNOptR/', purpose: 'Network optimization from data' },
            { name: 'GINsim', url: 'http://ginsim.org/', purpose: 'Gene regulatory network simulation' },
            { name: 'MaBoSS', url: 'https://maboss.curie.fr/', purpose: 'Boolean modeling of signaling' },
            { name: 'COPASI', url: 'http://copasi.org/', purpose: 'Biochemical system simulation' },
            { name: 'Virtual Cell', url: 'https://vcell.org/', purpose: 'Cell modeling platform' }
        ],
        reviews: [
            { title: 'Dynamics of cellular states', year: 2024, journal: 'Cell', doi: '10.1016/j.cell.2024.01.001' },
            { title: 'Cell signaling dynamics', year: 2023, journal: 'Nature Reviews MCB', doi: '10.1038/s41580-023-00599-x' }
        ]
    }
};