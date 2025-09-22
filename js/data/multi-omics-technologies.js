export const multiOmicsTechnologies = {
    id: 'multi-omics-technologies',
    name: 'Multi-Omics Technologies',
    description: 'Comprehensive molecular profiling technologies for understanding cellular states across genomics, transcriptomics, proteomics, and metabolomics layers',
    
    categories: [
        {
            id: 'genomics',
            name: 'Genomics & Epigenomics',
            icon: 'fas fa-dna',
            color: 'primary',
            items: [
                {
                    id: 'whole-genome-seq',
                    name: 'Single-Cell Whole Genome Sequencing',
                    description: 'Complete DNA sequence analysis at single-cell resolution to detect genetic variations, mutations, and copy number alterations',
                    specifications: {
                        coverage: '30-100x depth typical',
                        accuracy: '>99.9% base calling',
                        variants: 'SNVs, CNVs, SVs, indels',
                        resolution: 'Single nucleotide',
                        cellsPerRun: '100-1000 cells'
                    },
                    methods: [
                        'MALBAC (Multiple Annealing and Looping Based Amplification Cycles)',
                        'MDA (Multiple Displacement Amplification)',
                        'DOP-PCR (Degenerate Oligonucleotide PCR)',
                        'Primary template-directed amplification (PTA)',
                        'LIANTI (Linear Amplification via Transposon Insertion)'
                    ],
                    applications: [
                        'Cancer genomics - tracking clonal evolution',
                        'Prenatal genetic diagnosis from single cells',
                        'Studying mosaicism in neurodevelopment',
                        'Microbial single-cell genomics',
                        'Gamete genome analysis'
                    ],
                    keyMetrics: {
                        'Genome coverage': '85-95% typically achieved',
                        'Allelic dropout': '10-20% in best methods',
                        'False positive rate': '<1% for SNVs',
                        'Amplification bias': 'Reduced with newer methods like LIANTI'
                    },
                    databases: [
                        { name: 'gnomAD', url: 'https://gnomad.broadinstitute.org/', description: 'Genome Aggregation Database - population variant frequencies' },
                        { name: 'COSMIC', url: 'https://cancer.sanger.ac.uk/cosmic', description: 'Catalogue of Somatic Mutations in Cancer' },
                        { name: 'dbSNP', url: 'https://www.ncbi.nlm.nih.gov/snp/', description: 'Database of Single Nucleotide Polymorphisms' },
                        { name: 'ClinVar', url: 'https://www.ncbi.nlm.nih.gov/clinvar/', description: 'Clinical significance of genetic variants' }
                    ]
                },
                {
                    id: 'chromatin-accessibility',
                    name: 'Single-Cell Chromatin Accessibility',
                    description: 'Mapping open chromatin regions to understand gene regulatory landscapes and cell type-specific regulatory elements',
                    specifications: {
                        resolution: 'Single-cell with 1000-50000 peaks',
                        coverage: '1-10% of peaks detected per cell',
                        regulatory: 'Promoters, enhancers, silencers, insulators',
                        fragments: '1000-100000 per cell'
                    },
                    technologies: [
                        'scATAC-seq (10x Genomics Chromium)',
                        'sci-ATAC-seq (combinatorial indexing)',
                        'snATAC-seq (single nucleus ATAC)',
                        'scNOMe-seq (nucleosome occupancy + accessibility)',
                        'scMNase-seq (micrococcal nuclease)',
                        'scDNase-seq (DNase hypersensitivity)'
                    ],
                    challenges: [
                        'Extreme sparsity (~3% of sites detected)',
                        'Batch effects between experiments',
                        'Peak calling optimization needed',
                        'Doublet detection difficulties',
                        'Integration with gene expression data'
                    ],
                    analysis: ['ArchR', 'Signac', 'SnapATAC2', 'chromVAR', 'cisTopic', 'MAESTRO'],
                    biologicalInsights: [
                        'Cell type-specific regulatory elements',
                        'Transcription factor binding footprints',
                        'Chromatin remodeling dynamics',
                        'Enhancer-promoter interactions inference'
                    ],
                    databases: [
                        { name: 'ENCODE', url: 'https://www.encodeproject.org/', description: 'Encyclopedia of DNA Elements - chromatin accessibility maps' },
                        { name: 'Cistrome DB', url: 'http://cistrome.org/', description: 'Chromatin regulator ChIP-seq and chromatin accessibility data' }
                    ]
                },
                {
                    id: 'dna-methylation',
                    name: 'Single-Cell DNA Methylation',
                    description: 'Epigenetic modification patterns revealing cell identity, developmental history, and disease states',
                    coverage: '5-40% of CpGs per cell depending on method',
                    methods: [
                        'scBS-seq (single-cell bisulfite sequencing)',
                        'scRRBS (reduced representation bisulfite sequencing)',
                        'snmC-seq2 (single-nucleus methylcytosine sequencing)',
                        'sci-MET (combinatorial indexing methylation)',
                        'scNMT-seq (nucleosome, methylation, transcription)'
                    ],
                    biologicalRelevance: [
                        'Cell lineage tracing through methylation clocks',
                        'X-chromosome inactivation patterns',
                        'Imprinting control regions',
                        'Cancer methylation signatures',
                        'Aging and senescence markers'
                    ],
                    keyFacts: [
                        '~28 million CpGs in human genome',
                        '70-80% methylated in differentiated cells',
                        'CpG islands typically unmethylated',
                        'Methylation valleys mark developmental genes',
                        'PMDs (partially methylated domains) in cancer'
                    ],
                    databases: [
                        { name: 'MethBank', url: 'http://bigd.big.ac.cn/methbank/', description: 'DNA methylation data across species' },
                        { name: 'EWAS Data Hub', url: 'https://bigd.big.ac.cn/ewas/datahub', description: 'Epigenome-wide association studies' }
                    ]
                },
                {
                    id: 'histone-modifications',
                    name: 'Single-Cell Histone Modifications',
                    description: 'Chromatin state profiling through post-translational histone modifications',
                    marks: {
                        'H3K4me3': 'Active promoters, transcription initiation',
                        'H3K27ac': 'Active enhancers and promoters',
                        'H3K4me1': 'Enhancers (active and poised)',
                        'H3K36me3': 'Transcription elongation, gene bodies',
                        'H3K9me3': 'Constitutive heterochromatin',
                        'H3K27me3': 'Facultative heterochromatin, Polycomb silencing',
                        'H3K9ac': 'Active transcription'
                    },
                    techniques: [
                        'scCUT&Tag (Cleavage Under Targets and Tagmentation)',
                        'scCUT&RUN (Cleavage Under Targets and Release Using Nuclease)',
                        'scChIP-seq (chromatin immunoprecipitation)',
                        'DROP-ChIP (droplet microfluidics ChIP)',
                        'ACT-seq (antibody-guided chromatin tagmentation)',
                        'CoBATCH (combinatorial barcoding and targeted chromatin release)'
                    ],
                    resolution: '500-20000 peaks per cell',
                    applications: [
                        'Chromatin state mapping during differentiation',
                        'Enhancer identification and activity',
                        'Bivalent domain dynamics in stem cells',
                        'Heterochromatin formation in aging'
                    ],
                    databases: [
                        { name: 'ChIP-Atlas', url: 'https://chip-atlas.org/', description: 'Comprehensive ChIP-seq and epigenome database' },
                        { name: 'HistoneDB', url: 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4702785/', description: 'Histone sequences and variants' }
                    ]
                },
                {
                    id: '3d-genome',
                    name: 'Single-Cell 3D Genome Architecture',
                    description: 'Chromosome conformation and nuclear organization at single-cell resolution',
                    methods: [
                        'scHi-C (single-cell chromosome conformation capture)',
                        'sci-Hi-C (combinatorial indexing Hi-C)',
                        'scSPRITE (split-pool recognition of interactions)',
                        'scDam-C (DamID with Hi-C)',
                        'imaging methods (DNA FISH, ORCA)'
                    ],
                    structures: {
                        'TADs': 'Topologically associating domains (200kb-2Mb)',
                        'Loops': 'Enhancer-promoter interactions',
                        'Compartments': 'A (active) and B (inactive) compartments',
                        'LADs': 'Lamina-associated domains',
                        'Chromosomal territories': 'Spatial organization of chromosomes'
                    },
                    cellToCell: 'High variability in 3D structure between cells',
                    databases: [
                        { name: '4DN Data Portal', url: 'https://data.4dnucleome.org/', description: '4D Nucleome project - nuclear organization data' },
                        { name: '3DIV', url: 'http://3div.kr/', description: '3D-genome Interaction Viewer' }
                    ]
                }
            ]
        },
        {
            id: 'transcriptomics',
            name: 'Transcriptomics',
            icon: 'fas fa-scroll',
            color: 'secondary',
            items: [
                {
                    id: 'scrna-seq-comprehensive',
                    name: 'Single-Cell RNA Sequencing',
                    description: 'Comprehensive gene expression profiling revealing cell types, states, and dynamics',
                    platforms: {
                        'Droplet-based': {
                            '10x Chromium': 'Market leader, 500-10000 cells/channel',
                            'Drop-seq': 'Open-source, lower cost',
                            'inDrop': 'Hydrogel beads, CEL-Seq2 chemistry',
                            'ddSEQ': 'Bio-Rad system'
                        },
                        'Plate-based': {
                            'Smart-seq3': 'Full-length, high sensitivity, isoforms',
                            'MARS-seq': 'Massively parallel RNA single-cell',
                            'CEL-seq2': 'Early multiplexing method',
                            'Quartz-seq2': 'High sensitivity Japanese method'
                        },
                        'Combinatorial': {
                            'sci-RNA-seq': 'Split-pool barcoding, >1M cells',
                            'SPLiT-seq': 'No special equipment needed',
                            'SHARE-seq': 'RNA + chromatin simultaneously'
                        }
                    },
                    specifications: {
                        genes: '2000-6000 detected per cell typically',
                        cells: '100 to >1,000,000 per experiment',
                        sensitivity: '10-45% of transcripts captured',
                        accuracy: 'UMI-based removes PCR duplicates',
                        readDepth: '20,000-100,000 reads/cell'
                    },
                    biologicalInsights: [
                        'Novel cell type discovery',
                        'Continuous state transitions',
                        'Gene regulatory network inference',
                        'Cell-cell communication patterns',
                        'Disease mechanisms at cellular resolution'
                    ],
                    keyAdvances: [
                        'Cell hashing for sample multiplexing',
                        'Nuclear RNA sequencing for frozen tissues',
                        'Targeted gene panels (e.g., BD Rhapsody)',
                        'Direct RNA sequencing (no amplification)',
                        'In situ sequencing methods'
                    ],
                    databases: [
                        { name: 'Human Cell Atlas', url: 'https://www.humancellatlas.org/', description: 'Global effort to map all human cells' },
                        { name: 'Single Cell Portal', url: 'https://singlecell.broadinstitute.org/', description: 'Broad Institute single-cell data repository' },
                        { name: 'CellxGene', url: 'https://cellxgene.cziscience.com/', description: 'Chan Zuckerberg Initiative cell atlas' },
                        { name: 'PanglaoDB', url: 'https://panglaodb.se/', description: 'Single-cell RNA sequencing database' }
                    ]
                },
                {
                    id: 'spatial-transcriptomics',
                    name: 'Spatial Transcriptomics',
                    description: 'Gene expression with preserved spatial context for understanding tissue architecture',
                    categories: {
                        'Sequencing-based': [
                            'Visium (10x) - 55µm spots, whole transcriptome',
                            'Slide-seq V2 - 10µm beads, higher resolution',
                            'Stereo-seq - Adjustable 500nm-50µm resolution',
                            'DBiT-seq - Microfluidic barcoding, 10µm',
                            'Seq-Scope - Sub-micrometer resolution'
                        ],
                        'Imaging-based': [
                            'MERFISH - Error-robust, 100-10,000 genes',
                            'seqFISH+ - Up to 10,000 genes, true single-cell',
                            'ISS (In Situ Sequencing) - Padlock probes',
                            'ExSeq (Expansion Sequencing) - With tissue expansion',
                            'osmFISH - cyclic single molecule FISH'
                        ],
                        'Commercial platforms': [
                            'Xenium (10x) - In situ, subcellular, 5000 genes',
                            'CosMx SMI (NanoString) - 1000+ RNA, 64+ proteins',
                            'MERSCOPE (Vizgen) - MERFISH platform'
                        ]
                    },
                    resolution: {
                        'Visium': '55 µm spots, ~10 cells/spot',
                        'Visium HD': '2 µm bins, near single-cell',
                        'Slide-seq': '10 µm beads',
                        'MERFISH': '100-200 nm, subcellular',
                        'Stereo-seq': 'Flexible 500nm-50µm'
                    },
                    applications: [
                        'Tumor microenvironment mapping',
                        'Brain region molecular atlases',
                        'Developmental morphogen gradients',
                        'Disease pathology with spatial context',
                        'Cell-cell interaction networks'
                    ],
                    databases: [
                        { name: 'SpatialDB', url: 'https://www.spatialomics.org/SpatialDB/', description: 'Spatial transcriptomics database' },
                        { name: 'STOmics DB', url: 'https://db.cngb.org/stomics/', description: 'Spatio-temporal omics database' }
                    ]
                },
                {
                    id: 'rna-dynamics',
                    name: 'RNA Dynamics & Processing',
                    description: 'Temporal aspects of RNA metabolism revealing gene regulation dynamics',
                    methods: {
                        'RNA velocity': {
                            tools: 'scVelo, velocyto, VeloAE',
                            principle: 'Spliced vs unspliced RNA ratios',
                            timescale: 'Hours of future cell state'
                        },
                        'Metabolic labeling': {
                            techniques: 'sci-fate, NASC-seq, scSLAM-seq, TimeLapse-seq',
                            principle: '4sU or 6sG incorporation',
                            resolution: 'New vs old RNA, hours to days'
                        },
                        'RNA localization': {
                            methods: 'APEX-seq, FISH-seq, PrismEXT',
                            compartments: 'Nuclear, cytoplasmic, organellar'
                        },
                        'Alternative splicing': {
                            platforms: 'Smart-seq3, long-read sequencing',
                            tools: 'BRIE2, LeafCutter, Whippet',
                            detection: 'Isoform-level quantification'
                        }
                    },
                    biologicalProcesses: [
                        'Transcriptional bursting patterns',
                        'mRNA half-life determination',
                        'Translation efficiency measurement',
                        'RNA editing (A-to-I) detection',
                        'Circular RNA identification'
                    ],
                    timescales: {
                        'Immediate early genes': '15-30 minutes',
                        'Most mRNAs': '2-8 hour half-life',
                        'Long-lived mRNAs': 'Days (e.g., β-globin)',
                        'miRNAs': 'Days to weeks stability'
                    },
                    databases: [
                        { name: 'RADAR', url: 'http://rnaedit.com/', description: 'RNA editing database' },
                        { name: 'circBase', url: 'http://www.circbase.org/', description: 'Circular RNA database' }
                    ]
                },
                {
                    id: 'long-noncoding-rna',
                    name: 'Long Non-coding RNA',
                    description: 'Regulatory RNA molecules >200nt with diverse cellular functions',
                    characteristics: [
                        'Lower expression than mRNAs',
                        'Cell type-specific expression',
                        'Nuclear enrichment common',
                        'Poor sequence conservation'
                    ],
                    functions: {
                        'Chromatin regulation': 'XIST, HOTAIR, NEAT1',
                        'Transcriptional control': 'lincRNA-p21, MALAT1',
                        'Post-transcriptional': 'TINCR, BACE1-AS',
                        'Scaffolding': 'NEAT1 (paraspeckles), MALAT1'
                    },
                    detection: [
                        'Total RNA-seq required',
                        'Nuclear RNA enrichment helpful',
                        'Smart-seq better than droplet methods',
                        'CAGE-seq for 5\' ends'
                    ],
                    databases: [
                        { name: 'LNCipedia', url: 'https://lncipedia.org/', description: 'Comprehensive lncRNA database' },
                        { name: 'NONCODE', url: 'http://www.noncode.org/', description: 'Integrated non-coding RNA database' }
                    ]
                }
            ]
        },
        {
            id: 'proteomics',
            name: 'Proteomics & PTMs',
            icon: 'fas fa-cubes',
            color: 'success',
            items: [
                {
                    id: 'single-cell-proteomics',
                    name: 'Cell Proteomics',
                    description: 'Protein abundance and identification analysis',
                    technologies: {
                        'Mass spectrometry': {
                            'SCoPE-MS': 'Single Cell ProtEomics, 500-1000 proteins',
                            'nanoPOTS': 'Nanodroplet processing, improved sensitivity',
                            'OAD': 'One-pot analysis, simplified workflow',
                            'plexDIA': 'Multiplexed data-independent acquisition',
                            'Single-cell ProteomicsChip': 'Microfluidic integration'
                        },
                        'Antibody-based': {
                            'CyTOF': 'Mass cytometry, 40-50 markers',
                            'CITE-seq': 'RNA + 200+ surface proteins',
                            'REAP-seq': 'Similar to CITE-seq',
                            'AbSeq': 'BD Rhapsody protein detection',
                            'INs-seq': 'Intracellular protein + RNA'
                        },
                        'Imaging': {
                            'CODEX/PhenoCycler': '40-60 markers',
                            'MIBI-TOF': 'Ion beam imaging, 40+ markers',
                            '4i': 'Iterative immunofluorescence',
                            'CyCIF': 'Cyclic immunofluorescence',
                            'IBEX': 'Iterative bleaching extends multiplexing'
                        }
                    },
                    specifications: {
                        'MS-based': '500-3000 proteins/cell current',
                        'CyTOF': '40-50 simultaneous markers',
                        'CITE-seq': '200+ surface proteins + transcriptome',
                        'Imaging': '40-100 proteins with spatial info',
                        'Dynamic range': '7 orders of magnitude challenge'
                    },
                    challenges: [
                        'Limited sensitivity vs bulk methods',
                        'Protein copy number spans 10^7 range',
                        'No amplification possible',
                        'Post-translational modifications detection',
                        'Protein complex preservation'
                    ],
                    biologicalRelevance: [
                        'Direct functional readout',
                        'Signaling pathway states',
                        'Cell surface immunophenotype',
                        'Drug target expression',
                        'Post-transcriptional regulation'
                    ],
                    databases: [
                        { name: 'Human Protein Atlas', url: 'https://www.proteinatlas.org/', description: 'Protein expression across tissues and cells' },
                        { name: 'ProteomicsDB', url: 'https://www.proteomicsdb.org/', description: 'Multi-organism proteome resource' },
                        { name: 'PRIDE', url: 'https://www.ebi.ac.uk/pride/', description: 'PRoteomics IDEntifications database' },
                        { name: 'PeptideAtlas', url: 'http://www.peptideatlas.org/', description: 'Peptide and protein identification database' }
                    ]
                },
                {
                    id: 'ptm-profiling',
                    name: 'Post-Translational Modifications',
                    description: 'Chemical modifications that regulate protein function, localization, and stability',
                    types: {
                        'Phosphorylation': {
                            sites: '~230,000 known sites',
                            function: 'Signaling, activity regulation',
                            detection: 'Phospho-flow, MS with enrichment',
                            dynamics: 'Seconds to minutes'
                        },
                        'Ubiquitination': {
                            types: 'Mono, K48-linked, K63-linked',
                            function: 'Degradation, signaling, trafficking',
                            detection: 'Ubiquitin antibodies, MS',
                            complexity: '8 different linkage types'
                        },
                        'Acetylation': {
                            sites: '~20,000 known sites',
                            function: 'Gene regulation, metabolism',
                            enzymes: 'HATs and HDACs',
                            detection: 'Ac-lysine antibodies'
                        },
                        'Methylation': {
                            targets: 'Lysine, arginine residues',
                            states: 'Mono, di, tri-methylation',
                            function: 'Chromatin regulation, signaling',
                            writers: 'Methyltransferases'
                        },
                        'Glycosylation': {
                            types: 'N-glycosylation, O-glycosylation',
                            prevalence: '>50% of proteins',
                            function: 'Folding, recognition, stability',
                            complexity: 'Enormous structural diversity'
                        }
                    },
                    methods: [
                        'Phospho-specific flow cytometry',
                        'PTM-specific antibodies in CyTOF',
                        'MS with enrichment strategies',
                        'Proximity ligation assays',
                        'FRET-based biosensors'
                    ],
                    dynamics: {
                        'Phosphorylation': 'Seconds to minutes',
                        'Ubiquitination': 'Minutes to hours',
                        'Acetylation': 'Minutes to hours',
                        'Methylation': 'Hours to days',
                        'Glycosylation': 'Hours (ER/Golgi processing)'
                    },
                    databases: [
                        { name: 'PhosphoSitePlus', url: 'https://www.phosphosite.org/', description: 'PTMs, especially phosphorylation' },
                        { name: 'PTMcode', url: 'https://ptmcode.embl.de/', description: 'Functional associations of PTMs' },
                        { name: 'GlyGen', url: 'https://www.glygen.org/', description: 'Glycosylation data integration' },
                        { name: 'UniMod', url: 'http://www.unimod.org/', description: 'Protein modifications for mass spectrometry' }
                    ]
                },
                {
                    id: 'protein-interactions',
                    name: 'Protein-Protein Interactions',
                    description: 'Mapping molecular interaction networks that drive cellular functions',
                    methods: [
                        'Proximity labeling (BioID, APEX, TurboID)',
                        'Split-protein reconstitution (BiFC, split-GFP)',
                        'FRET/BRET for dynamics',
                        'Co-immunoprecipitation MS',
                        'Cross-linking MS (XL-MS)',
                        'Protein complementation assays'
                    ],
                    scale: {
                        'Human interactome': '~650,000 possible interactions',
                        'Core interactions': '~15,000-25,000 high-confidence',
                        'Dynamic range': 'pM to mM binding affinities',
                        'Proteins per complex': '2 to >100 subunits'
                    },
                    types: [
                        'Stable complexes (ribosomes, proteasomes)',
                        'Transient interactions (signaling)',
                        'Enzyme-substrate interactions',
                        'Scaffold-mediated assemblies',
                        'Phase-separated condensates'
                    ],
                    singleCellApproaches: [
                        'PLA (Proximity Ligation Assay)',
                        'Single-molecule FRET',
                        'Multiplexed imaging of complexes',
                        'Flow cytometry with interaction readouts'
                    ],
                    databases: [
                        { name: 'STRING', url: 'https://string-db.org/', description: 'Protein-protein interaction networks' },
                        { name: 'BioGRID', url: 'https://thebiogrid.org/', description: 'Genetic and protein interactions' },
                        { name: 'IntAct', url: 'https://www.ebi.ac.uk/intact/', description: 'Molecular interaction database' },
                        { name: 'CORUM', url: 'http://mips.helmholtz-muenchen.de/corum/', description: 'Comprehensive resource of mammalian protein complexes' }
                    ]
                },
                {
                    id: 'protein-localization',
                    name: 'Protein Localization',
                    description: 'Subcellular protein distribution and trafficking',
                    compartments: [
                        'Nuclear vs cytoplasmic',
                        'Membrane-bound organelles',
                        'Membrane proteins (PM, ER, Golgi)',
                        'Cytoskeletal association',
                        'Phase-separated condensates'
                    ],
                    methods: [
                        'Immunofluorescence microscopy',
                        'Protein-fragment complementation',
                        'APEX/HRP proximity labeling',
                        'Fractionation + MS',
                        'Live-cell imaging with FPs'
                    ],
                    dynamicTrafficking: [
                        'ER to Golgi transport',
                        'Nuclear import/export',
                        'Endocytosis and recycling',
                        'Mitochondrial import',
                        'Peroxisomal targeting'
                    ],
                    databases: [
                        { name: 'Human Protein Atlas - Subcellular', url: 'https://www.proteinatlas.org/subcellular', description: 'Subcellular protein localization' },
                        { name: 'COMPARTMENTS', url: 'https://compartments.jensenlab.org/', description: 'Subcellular localization database' }
                    ]
                }
            ]
        },
        {
            id: 'metabolomics',
            name: 'Metabolomics & Lipidomics',
            icon: 'fas fa-flask',
            color: 'warning',
            items: [
                {
                    id: 'single-cell-metabolomics',
                    name: 'Single-Cell Metabolomics',
                    description: 'Small molecule profiling revealing metabolic states and activities in individual cells',
                    technologies: [
                        'Live single-cell MS (Live-seq MS)',
                        'MALDI-MS imaging (spatial metabolomics)',
                        'Secondary ion MS (SIMS) - subcellular resolution',
                        'Raman spectroscopy - label-free',
                        'Capillary electrophoresis MS',
                        'Acoustic droplet ejection MS'
                    ],
                    metabolites: {
                        'Primary': 'Amino acids, nucleotides, sugars, organic acids',
                        'Secondary': 'Lipids, cofactors, vitamins, hormones',
                        'Coverage': '50-500 metabolites detectable per cell',
                        'Concentration': 'pM to mM range in cells'
                    },
                    challenges: [
                        'Rapid metabolite turnover (seconds)',
                        'Chemical diversity requiring multiple methods',
                        'Matrix effects and ion suppression',
                        'Sample preparation artifacts',
                        'Limited amplification possibilities'
                    ],
                    applications: [
                        'Cancer metabolism heterogeneity',
                        'Drug metabolism and resistance',
                        'Metabolic disease mechanisms',
                        'Immune cell activation states',
                        'Stem cell metabolism'
                    ],
                    databases: [
                        { name: 'HMDB', url: 'https://hmdb.ca/', description: 'Human Metabolome Database' },
                        { name: 'METLIN', url: 'https://metlin.scripps.edu/', description: 'Metabolite and tandem MS database' },
                        { name: 'MetaboLights', url: 'https://www.ebi.ac.uk/metabolights/', description: 'Metabolomics experiments and metadata' }
                    ]
                },
                {
                    id: 'metabolic-flux',
                    name: 'Metabolic Flux Analysis',
                    description: 'Dynamic measurement of metabolite flow through biochemical pathways',
                    methods: [
                        '13C isotope tracing (glucose, glutamine)',
                        '2H (deuterium) labeling',
                        'Seahorse analysis (OCR/ECAR)',
                        'FLIM-FRET metabolic sensors',
                        'Two-photon FLIM (NAD(P)H)',
                        'Single-cell isotope tracing'
                    ],
                    pathways: [
                        'Glycolysis (Warburg effect in cancer)',
                        'TCA cycle (oxidative metabolism)',
                        'Pentose phosphate pathway (NADPH, ribose)',
                        'Fatty acid synthesis/oxidation',
                        'One-carbon metabolism (methylation, nucleotides)',
                        'Amino acid metabolism'
                    ],
                    timescales: {
                        'Glycolysis': 'Seconds to minutes',
                        'TCA cycle': 'Minutes',
                        'Lipogenesis': 'Hours',
                        'Protein synthesis': 'Hours to days'
                    },
                    cellularReadouts: [
                        'ATP/ADP/AMP ratios',
                        'NAD+/NADH redox state',
                        'NADP+/NADPH for biosynthesis',
                        'ROS levels and oxidative stress',
                        'pH changes',
                        'Mitochondrial membrane potential'
                    ],
                    databases: [
                        { name: 'KEGG', url: 'https://www.kegg.jp/', description: 'Metabolic pathway maps' },
                        { name: 'Reactome', url: 'https://reactome.org/', description: 'Pathway database including metabolism' }
                    ]
                },
                {
                    id: 'lipidomics',
                    name: 'Lipidomics',
                    description: 'Comprehensive analysis of cellular lipids and their biological roles',
                    classes: {
                        'Structural': {
                            types: 'Phospholipids, sphingolipids, cholesterol',
                            function: 'Membrane structure and fluidity'
                        },
                        'Storage': {
                            types: 'Triglycerides, cholesterol esters',
                            function: 'Energy storage, droplet formation'
                        },
                        'Signaling': {
                            types: 'Eicosanoids, phosphoinositides, ceramides',
                            function: 'Cell signaling and regulation'
                        },
                        'Complex': {
                            types: 'Cardiolipins, gangliosides',
                            function: 'Organelle-specific functions'
                        }
                    },
                    diversity: '>40,000 unique lipid species possible',
                    methods: [
                        'LC-MS/MS (liquid chromatography MS)',
                        'Shotgun lipidomics (direct infusion)',
                        'Ion mobility spectrometry (IMS)',
                        'MALDI-MS imaging',
                        'Desorption ESI (DESI)',
                        'Single-cell lipidomics'
                    ],
                    cellularRoles: [
                        'Membrane composition and domains',
                        'Energy storage and mobilization',
                        'Signaling cascade initiation',
                        'Protein anchoring and trafficking',
                        'Inflammation and immunity'
                    ],
                    databases: [
                        { name: 'LIPID MAPS', url: 'https://www.lipidmaps.org/', description: 'Comprehensive lipid database' },
                        { name: 'SwissLipids', url: 'https://www.swisslipids.org/', description: 'Curated lipid knowledge resource' },
                        { name: 'LipidBank', url: 'http://lipidbank.jp/', description: 'Japanese lipid database' }
                    ]
                },
                {
                    id: 'immunometabolism',
                    name: 'Immunometabolism',
                    description: 'Metabolic control of immune cell function and fate',
                    states: {
                        'Naive T cells': 'OXPHOS, low metabolic rate',
                        'Activated T cells': 'Aerobic glycolysis, high anabolism',
                        'Memory T cells': 'FAO, OXPHOS, self-renewal',
                        'M1 macrophages': 'Glycolysis, pro-inflammatory',
                        'M2 macrophages': 'OXPHOS, anti-inflammatory'
                    },
                    metabolicCheckpoints: [
                        'mTOR signaling in T cell activation',
                        'AMPK in metabolic stress',
                        'HIF-1α in hypoxic response',
                        'Metabolite regulation of epigenetics'
                    ],
                    therapeuticTargets: [
                        'Glycolysis inhibitors in autoimmunity',
                        'IDO inhibitors in cancer',
                        'Glutamine metabolism in tumors',
                        'Itaconate as anti-inflammatory'
                    ]
                }
            ]
        },
        {
            id: 'multimodal-integration',
            name: 'Multimodal Integration',
            icon: 'fas fa-layer-group',
            color: 'info',
            items: [
                {
                    id: 'paired-omics',
                    name: 'Paired Multi-Omics Technologies',
                    description: 'Simultaneous measurement of multiple molecular modalities from the same single cells',
                    technologies: {
                        'CITE-seq': {
                            modalities: 'RNA + surface proteins',
                            scale: 'Thousands of cells, 200+ proteins',
                            vendor: 'New York Genome Center method'
                        },
                        '10x Multiome': {
                            modalities: 'RNA + chromatin accessibility',
                            scale: 'Thousands of cells simultaneously',
                            vendor: '10x Genomics'
                        },
                        'SHARE-seq': {
                            modalities: 'RNA + chromatin accessibility',
                            scale: 'Hundreds of thousands of cells',
                            vendor: 'Academic method'
                        },
                        'DOGMA-seq': {
                            modalities: 'RNA + ATAC + proteins',
                            scale: 'Triple modality',
                            vendor: 'NYGC/Memorial Sloan Kettering'
                        },
                        'TEA-seq': {
                            modalities: 'RNA + ATAC + proteins + TCR/BCR',
                            scale: 'Quad modality for immune cells',
                            vendor: 'Academic development'
                        },
                        'NEAT-seq': {
                            modalities: 'RNA + proteins + accessible chromatin',
                            scale: 'Nuclear and cytoplasmic',
                            vendor: 'Recent development'
                        }
                    },
                    advantages: [
                        'True single-cell correlations',
                        'Reduced batch effects',
                        'Direct regulatory relationships',
                        'Comprehensive cell state definition',
                        'Lineage and state coupling'
                    ],
                    computationalMethods: [
                        'WNN (Weighted Nearest Neighbor) - Seurat',
                        'MOFA+ (Multi-Omics Factor Analysis)',
                        'totalVI (total Variational Inference)',
                        'Cobolt (Bayesian integration)',
                        'MultiVI (deep generative modeling)',
                        'GLUE (Graph-linked unified embedding)'
                    ],
                    databases: [
                        { name: 'Single Cell Portal', url: 'https://singlecell.broadinstitute.org/', description: 'Multi-modal single-cell data' },
                        { name: 'HuBMAP', url: 'https://hubmapconsortium.org/', description: 'Human BioMolecular Atlas Program' }
                    ]
                },
                // {
                //     id: 'paired-omics',
                //     name: 'Paired Multi-Omics Technologies',
                //     description: 'Simultaneous measurement of multiple molecular modalities from the same single cells',
                //     technologies: [
                //         {
                //             name: 'CITE-seq',
                //             modalities: 'RNA + surface proteins',
                //             scale: 'Thousands of cells, 200+ proteins',
                //             vendor: 'New York Genome Center method'
                //         },
                //         {
                //             name: '10x Multiome',
                //             modalities: 'RNA + chromatin accessibility',
                //             scale: 'Thousands of cells simultaneously',
                //             vendor: '10x Genomics'
                //         },
                //         {
                //             name: 'SHARE-seq',
                //             modalities: 'RNA + chromatin accessibility',
                //             scale: 'Hundreds of thousands of cells',
                //             vendor: 'Academic method'
                //         },
                //         {
                //             name: 'DOGMA-seq',
                //             modalities: 'RNA + ATAC + proteins',
                //             scale: 'Triple modality',
                //             vendor: 'NYGC/Memorial Sloan Kettering'
                //         },
                //         {
                //             name: 'TEA-seq',
                //             modalities: 'RNA + ATAC + proteins + TCR/BCR',
                //             scale: 'Quad modality for immune cells',
                //             vendor: 'Academic development'
                //         },
                //         {
                //             name: 'NEAT-seq',
                //             modalities: 'RNA + proteins + accessible chromatin',
                //             scale: 'Nuclear and cytoplasmic',
                //             vendor: 'Recent development'
                //         }
                //     ],
                //     advantages: [
                //         'True single-cell correlations',
                //         'Reduced batch effects',
                //         'Direct regulatory relationships',
                //         'Comprehensive cell state definition',
                //         'Lineage and state coupling'
                //     ],
                //     computationalMethods: [
                //         'WNN (Weighted Nearest Neighbor) - Seurat',
                //         'MOFA+ (Multi-Omics Factor Analysis)',
                //         'totalVI (total Variational Inference)',
                //         'Cobolt (Bayesian integration)',
                //         'MultiVI (deep generative modeling)',
                //         'GLUE (Graph-linked unified embedding)'
                //     ],
                //     databases: [
                //         { name: 'Single Cell Portal', url: 'https://singlecell.broadinstitute.org/', description: 'Multi-modal single-cell data' },
                //         { name: 'HuBMAP', url: 'https://hubmapconsortium.org/', description: 'Human BioMolecular Atlas Program' }
                //     ]
                // },
                {
                    id: 'temporal-profiling',
                    name: 'Temporal Multi-Omics',
                    description: 'Time-resolved molecular dynamics across multiple layers',
                    approaches: [
                        'Metabolic labeling (4sU-seq, SLAM-seq, TUC-seq)',
                        'Cell cycle synchronization and staging',
                        'Lineage tracing (CRISPR barcoding, lentiviral)',
                        'Live imaging followed by omics',
                        'Pulse-chase experiments',
                        'Developmental time series'
                    ],
                    timescales: {
                        'Phosphorylation': 'Seconds to minutes',
                        'Transcription': 'Minutes to hours',
                        'Translation': 'Hours',
                        'Protein degradation': 'Hours to days',
                        'Epigenetic changes': 'Hours to generations'
                    },
                    applications: [
                        'Cell differentiation trajectories',
                        'Drug response kinetics',
                        'Circadian rhythm studies',
                        'Disease progression modeling',
                        'Developmental dynamics',
                        'Aging processes'
                    ],
                    challenges: [
                        'Synchronization of processes',
                        'Destructive measurement',
                        'Computational trajectory inference',
                        'Batch effect over time'
                    ],
                    tools: [
                        'Monocle3 (trajectory inference)',
                        'scVelo (RNA velocity)',
                        'CellRank (fate mapping)',
                        'Tempora (temporal ordering)'
                    ]
                },
                {
                    id: 'spatial-multiomics',
                    name: 'Spatial Multi-Omics',
                    description: 'Multiple molecular modalities with spatial information',
                    technologies: [
                        'DBiT-seq (RNA + proteins spatial)',
                        'Spatial-CITE-seq',
                        'SM-Omics (RNA + proteins)',
                        'Spatial epigenome-transcriptome co-profiling',
                        'MERFISH + protein imaging'
                    ],
                    integration: [
                        'Registration across modalities',
                        'Cell segmentation challenges',
                        'Multi-scale integration',
                        'Neighborhood analysis'
                    ],
                    biologicalInsights: [
                        'Tissue organization principles',
                        'Cell-cell communication in space',
                        'Microenvironment influence',
                        'Morphogen gradient responses'
                    ],
                    databases: [
                        { name: 'Spatial Omics Database', url: 'https://www.spatialomics.org/', description: 'Spatial multi-omics repository' }
                    ]
                }
            ]
        }
    ],
    
    resources: {
        databases: [
            { name: 'Human Cell Atlas', url: 'https://www.humancellatlas.org/', type: 'Comprehensive single-cell reference' },
            { name: 'ENCODE', url: 'https://www.encodeproject.org/', type: 'Encyclopedia of DNA Elements' },
            { name: 'GTEx', url: 'https://gtexportal.org/', type: 'Genotype-Tissue Expression' },
            { name: 'HuBMAP', url: 'https://hubmapconsortium.org/', type: 'Human BioMolecular Atlas' },
            { name: 'CZ CELLxGENE', url: 'https://cellxgene.cziscience.com/', type: 'Single-cell data platform' },
            { name: '4DN', url: 'https://data.4dnucleome.org/', type: '4D Nucleome - nuclear organization' }
        ],
        tools: [
            { name: 'Seurat', url: 'https://satijalab.org/seurat/', purpose: 'R toolkit for single-cell analysis' },
            { name: 'Scanpy', url: 'https://scanpy.readthedocs.io/', purpose: 'Python single-cell analysis' },
            { name: 'MOFA+', url: 'https://biofam.github.io/MOFA2/', purpose: 'Multi-omics factor analysis' },
            { name: 'mixOmics', url: 'http://mixomics.org/', purpose: 'Multivariate omics integration' }
        ],
        // reviews: [
        //     { title: 'Integrative single-cell analysis', year: 2024, journal: 'Nature Reviews Genetics', doi: '10.1038/s41576-024-00716-y' },
        //     { title: 'Single-cell multi-omics', year: 2023, journal: 'Nature Methods', doi: '10.1038/s41592-023-01992-y' }
        // ]
        reviews: [
            {
                title: "Best practices for single-cell analysis across modalities",
                year: 2023,
                journal: "Nature Reviews Genetics",
                doi: "10.1038/s41576-023-00586-w"
            },
            {
                title: "Advances in single-cell omics and multiomics for high-resolution molecular profiling",
                year: 2024,
                journal: "Experimental & Molecular Medicine",
                doi: "10.1038/s12276-024-01186-2"
            },
            {
                title: "Spatial Transcriptomics Brings New Challenges and Opportunities",
                year: 2024,
                journal: "Annual Review of Biomedical Data Science",
                doi: "10.1146/annurev-biodatasci-040324-030052"
            },
            {
                title: "Spatial and single-cell profiling of the metabolome",
                year: 2023,
                journal: "Nature Metabolism",
                doi: "10.1038/s43587-023-00513-y"
            },
            {
                title: "Single-cell proteomics enabled by next-generation sequencing or mass spectrometry",
                year: 2023,
                journal: "Nature Methods",
                doi: "10.1038/s41592-023-01791-5"
            },
            {
                title: "Systematic benchmarking of single-cell ATAC-sequencing protocols",
                year: 2024,
                journal: "Nature Biotechnology",
                doi: "10.1038/s41587-023-01881-x"
            },
            {
                title: "Long non-coding RNAs: definitions, functions, challenges and recommendations",
                year: 2023,
                journal: "Nature Reviews Molecular Cell Biology",
                doi: "10.1038/s41580-022-00566-8"
            },
            {
                title: "Benchmarking multi-omics integration algorithms across multiple omics types",
                year: 2024,
                journal: "Briefings in Bioinformatics",
                doi: "10.1093/bib/bbae095"
            },
            {
                title: "Integration of spatial and single-cell data across modalities with MaxFuse",
                year: 2024,
                journal: "Nature Biotechnology",
                doi: "10.1038/s41587-023-01935-0"
            },
            {
                title: "SpatialGlue deciphers spatial domains from spatial multi-omics",
                year: 2024,
                journal: "Nature Methods",
                doi: "10.1038/s41592-024-02316-4"
            }
        ]
    }
};