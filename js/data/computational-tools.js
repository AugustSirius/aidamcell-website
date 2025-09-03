export const computationalTools = {
    id: 'computational-tools',
    name: 'Computational Tools',
    description: 'Software, algorithms, and analysis methods for single-cell data processing, analysis, and interpretation',
    
    categories: [
        {
            id: 'preprocessing',
            name: 'Preprocessing & QC',
            icon: 'fas fa-filter',
            color: 'primary',
            items: [
                {
                    id: 'cellranger',
                    name: 'Cell Ranger',
                    description: '10x Genomics official pipeline for processing Chromium single-cell data',
                    specifications: {
                        language: 'Binary (Python/Rust)',
                        license: 'MIT (v7.0+), Commercial (earlier)',
                        platform: 'Linux/Mac',
                        inputFormat: 'BCL/FASTQ',
                        outputFormat: 'Feature-barcode matrix, BAM, cloupe'
                    },
                    features: [
                        'Barcode processing and correction',
                        'Alignment with STAR',
                        'UMI counting and deduplication',
                        'Cell calling with EmptyDrops',
                        'Feature barcoding support',
                        'Intron mode for nuclei',
                        'Secondary analysis (PCA, clustering, t-SNE)'
                    ],
                    workflow: [
                        'mkfastq - Convert BCL to FASTQ',
                        'count - Align reads and generate feature-barcode matrices',
                        'aggr - Aggregate multiple samples with batch correction',
                        'reanalyze - Rerun secondary analysis with custom parameters',
                        'multi - Process multiplexed samples with hashtags'
                    ],
                    advantages: [
                        'Optimized for 10x Genomics data',
                        'Comprehensive documentation',
                        'Regular updates and support',
                        'Integrated Loupe Browser visualization'
                    ],
                    limitations: [
                        'Primarily designed for 10x data',
                        'Limited customization options',
                        'High memory requirements (64GB+ recommended)',
                        'Closed-source components'
                    ],
                    databases: [
                        { name: '10x Support', url: 'https://support.10xgenomics.com/', description: 'Official documentation and tutorials' },
                        { name: '10x Downloads', url: 'https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest', description: 'Software downloads' }
                    ]
                },
                {
                    id: 'starsolo',
                    name: 'STARsolo',
                    description: 'STAR aligner with built-in single-cell RNA-seq capabilities',
                    specifications: {
                        language: 'C++',
                        license: 'MIT',
                        platform: 'Linux/Mac/Windows',
                        speed: '2-5x faster than Cell Ranger'
                    },
                    features: [
                        'Multiple barcode geometries support',
                        'Cell calling algorithms (EmptyDrops, simple)',
                        'Gene expression and velocity output',
                        'Multimodal data processing',
                        'Custom barcode whitelist'
                    ],
                    advantages: [
                        'Highly customizable',
                        'Fast processing speed',
                        'Supports various protocols',
                        'Active development'
                    ],
                    databases: [
                        { name: 'STAR GitHub', url: 'https://github.com/alexdobin/STAR', description: 'Source code and documentation' },
                        { name: 'STAR Manual', url: 'https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf', description: 'Comprehensive manual' }
                    ]
                },
                {
                    id: 'kallisto-bustools',
                    name: 'kallisto | bustools',
                    description: 'Pseudoalignment-based fast single-cell RNA-seq quantification',
                    specifications: {
                        language: 'C++/Python',
                        license: 'BSD-2',
                        speed: '10-100x faster than alignment',
                        accuracy: 'Comparable to alignment methods'
                    },
                    features: [
                        'Near-optimal pseudoalignment',
                        'Barcode error correction',
                        'UMI collapsing',
                        'RNA velocity preprocessing',
                        'Multiple technology support'
                    ],
                    workflow: [
                        'kallisto bus - Generate BUS files',
                        'bustools correct - Error correction',
                        'bustools sort - Sort BUS files',
                        'bustools count - Generate count matrix'
                    ],
                    databases: [
                        { name: 'kallistobus.tools', url: 'https://www.kallistobus.tools/', description: 'Official tutorials and workflows' }
                    ]
                },
                {
                    id: 'alevin-fry',
                    name: 'alevin-fry',
                    description: 'Fast and accurate single-cell quantification using selective alignment',
                    specifications: {
                        language: 'Rust/C++',
                        license: 'BSD-3',
                        method: 'Selective alignment',
                        memory: 'Low memory footprint'
                    },
                    features: [
                        'Selective alignment strategy',
                        'Advanced UMI resolution',
                        'Barcode correction',
                        'Spliced/unspliced quantification',
                        'Sketch-based indexing'
                    ],
                    databases: [
                        { name: 'alevin-fry docs', url: 'https://alevin-fry.readthedocs.io/', description: 'Documentation and tutorials' }
                    ]
                },
                {
                    id: 'scrublet',
                    name: 'Scrublet',
                    description: 'Computational doublet detection in single-cell transcriptomic data',
                    specifications: {
                        language: 'Python',
                        method: 'Simulation-based',
                        accuracy: '85-95% depending on dataset',
                        FDR: '~5-10%'
                    },
                    features: [
                        'Simulated doublet generation',
                        'K-nearest neighbor graph',
                        'Doublet score calculation',
                        'Automatic threshold detection'
                    ],
                    databases: [
                        { name: 'Scrublet GitHub', url: 'https://github.com/AllonKleinLab/scrublet', description: 'Source code and examples' }
                    ]
                },
                {
                    id: 'doubletfinder',
                    name: 'DoubletFinder',
                    description: 'R-based doublet detection optimized for different scRNA-seq protocols',
                    specifications: {
                        language: 'R',
                        integration: 'Seurat workflow',
                        method: 'pN-pK parameter sweep',
                        accuracy: '>90% for 10x data'
                    },
                    features: [
                        'Protocol-specific optimization',
                        'Parameter sweep automation',
                        'Ground-truth informed',
                        'Visualization tools'
                    ],
                    databases: [
                        { name: 'DoubletFinder GitHub', url: 'https://github.com/chris-mcginnis-ucsf/DoubletFinder', description: 'Implementation and tutorials' }
                    ]
                },
                {
                    id: 'soupx',
                    name: 'SoupX',
                    description: 'Ambient RNA removal from droplet-based single-cell RNA-seq data',
                    specifications: {
                        language: 'R',
                        method: 'Contamination fraction estimation',
                        compatibility: '10x, Drop-seq, inDrop'
                    },
                    features: [
                        'Automatic contamination estimation',
                        'Manual marker genes option',
                        'Visual QC plots',
                        'Count matrix correction'
                    ],
                    databases: [
                        { name: 'SoupX GitHub', url: 'https://github.com/constantAmateur/SoupX', description: 'Package and documentation' }
                    ]
                }
            ]
        },
        {
            id: 'analysis-frameworks',
            name: 'Analysis Frameworks',
            icon: 'fas fa-chart-line',
            color: 'secondary',
            items: [
                {
                    id: 'seurat',
                    name: 'Seurat',
                    description: 'Comprehensive R toolkit for single-cell genomics analysis',
                    specifications: {
                        language: 'R',
                        version: '5.0+ (latest)',
                        license: 'MIT',
                        citation: 'Hao et al., Cell 2021; Stuart et al., Cell 2019'
                    },
                    features: [
                        'Quality control and filtering',
                        'Normalization (LogNormalize, SCTransform v2)',
                        'Data integration (CCA, RPCA, Harmony, FastMNN)',
                        'Dimension reduction (PCA, ICA, NMF)',
                        'Clustering (Louvain, Leiden, SLM)',
                        'Differential expression (multiple methods)',
                        'Trajectory inference integration',
                        'Spatial data analysis',
                        'Multimodal integration (WNN)',
                        'Cell type annotation'
                    ],
                    workflow: [
                        'CreateSeuratObject - Initialize data structure',
                        'PercentageFeatureSet - Calculate QC metrics',
                        'NormalizeData/SCTransform - Normalization',
                        'FindVariableFeatures - Feature selection',
                        'ScaleData - Scaling and centering',
                        'RunPCA - Linear dimension reduction',
                        'FindNeighbors - KNN graph construction',
                        'FindClusters - Graph-based clustering',
                        'RunUMAP/RunTSNE - Non-linear embedding',
                        'FindAllMarkers - Marker identification'
                    ],
                    advantages: [
                        'Most comprehensive feature set',
                        'Excellent documentation and vignettes',
                        'Large active community',
                        'Regular updates and maintenance',
                        'Extensive visualization options',
                        'Good memory management with v5'
                    ],
                    limitations: [
                        'R memory constraints for very large datasets',
                        'Steep learning curve for beginners',
                        'Some operations can be slow'
                    ],
                    databases: [
                        { name: 'Seurat Website', url: 'https://satijalab.org/seurat/', description: 'Official tutorials and vignettes' },
                        { name: 'Seurat GitHub', url: 'https://github.com/satijalab/seurat', description: 'Source code and issues' },
                        { name: 'Seurat Data', url: 'https://github.com/satijalab/seurat-data', description: 'Example datasets' }
                    ]
                },
                {
                    id: 'scanpy',
                    name: 'Scanpy',
                    description: 'Scalable Python toolkit for single-cell analysis',
                    specifications: {
                        language: 'Python',
                        license: 'BSD-3',
                        dataStructure: 'AnnData',
                        gpu: 'RAPIDS support available',
                        citation: 'Wolf et al., Genome Biology 2018'
                    },
                    features: [
                        'Preprocessing and quality control',
                        'Highly variable gene selection',
                        'Batch correction (Combat, BBKNN, Scanorama)',
                        'Clustering (Louvain, Leiden)',
                        'Trajectory inference (PAGA, DPT)',
                        'Differential expression',
                        'Gene set enrichment',
                        'Integration with scverse ecosystem'
                    ],
                    ecosystem: [
                        'AnnData - Annotated data matrix',
                        'Squidpy - Spatial analysis',
                        'scVelo - RNA velocity',
                        'CellRank - Fate mapping',
                        'Pertpy - Perturbation analysis',
                        'muon - Multimodal analysis'
                    ],
                    advantages: [
                        'Python ecosystem integration',
                        'Memory efficient with backed mode',
                        'GPU acceleration support',
                        'Modular and extensible',
                        'Excellent plotting capabilities',
                        'Active development'
                    ],
                    databases: [
                        { name: 'Scanpy Docs', url: 'https://scanpy.readthedocs.io/', description: 'Documentation and tutorials' },
                        { name: 'scverse', url: 'https://scverse.org/', description: 'Single-cell analysis ecosystem' }
                    ]
                },
                {
                    id: 'bioconductor-suite',
                    name: 'Bioconductor Single-Cell',
                    description: 'Collection of R packages for rigorous single-cell analysis',
                    specifications: {
                        language: 'R',
                        license: 'Various open source',
                        framework: 'SingleCellExperiment',
                        packages: '100+ packages'
                    },
                    keyPackages: [
                        'SingleCellExperiment - Data infrastructure',
                        'scater - Visualization and QC',
                        'scran - Normalization and batch correction',
                        'DropletUtils - Droplet data processing',
                        'batchelor - Batch effect removal',
                        'bluster - Clustering algorithms',
                        'scDblFinder - Doublet detection'
                    ],
                    advantages: [
                        'Rigorous statistical methods',
                        'Excellent documentation standards',
                        'Interoperable ecosystem',
                        'Peer-reviewed packages',
                        'Strong developer community'
                    ],
                    databases: [
                        { name: 'Bioconductor', url: 'https://bioconductor.org/packages/devel/workflows/vignettes/simpleSingleCell/inst/doc/', description: 'Single-cell workflow' },
                        { name: 'OSCA', url: 'https://bioconductor.org/books/release/OSCA/', description: 'Orchestrating Single-Cell Analysis book' }
                    ]
                },
                {
                    id: 'scvi-tools',
                    name: 'scvi-tools',
                    description: 'Deep learning methods for single-cell omics',
                    specifications: {
                        language: 'Python',
                        framework: 'PyTorch',
                        license: 'BSD-3',
                        gpu: 'CUDA acceleration'
                    },
                    models: [
                        'scVI - Variational inference for scRNA-seq',
                        'scANVI - Semi-supervised annotation',
                        'totalVI - CITE-seq analysis',
                        'PeakVI - scATAC-seq analysis',
                        'MultiVI - Multimodal integration',
                        'DestVI - Spatial deconvolution'
                    ],
                    advantages: [
                        'Principled probabilistic models',
                        'Uncertainty quantification',
                        'Scalable to millions of cells',
                        'Transfer learning capabilities',
                        'GPU acceleration'
                    ],
                    databases: [
                        { name: 'scvi-tools', url: 'https://scvi-tools.org/', description: 'Documentation and tutorials' }
                    ]
                }
            ]
        },
        {
            id: 'specialized-tools',
            name: 'Specialized Analysis',
            icon: 'fas fa-tools',
            color: 'success',
            items: [
                {
                    id: 'monocle3',
                    name: 'Monocle 3',
                    description: 'Trajectory inference and pseudotime analysis toolkit',
                    specifications: {
                        language: 'R',
                        method: 'UMAP + principal graph (reversed graph embedding)',
                        citation: 'Cao et al., Nature 2019',
                        scalability: 'Millions of cells'
                    },
                    features: [
                        'Trajectory inference without prior knowledge',
                        'Pseudotime ordering',
                        'Branch point detection and analysis',
                        'Differential expression along trajectories',
                        'Gene module discovery',
                        'Subclustering along paths'
                    ],
                    workflow: [
                        'new_cell_data_set - Create CDS object',
                        'preprocess_cds - Normalization and PCA',
                        'reduce_dimension - UMAP dimension reduction',
                        'cluster_cells - Optional clustering',
                        'learn_graph - Learn principal graph',
                        'order_cells - Assign pseudotime',
                        'graph_test - Test for trajectory-dependent expression'
                    ],
                    advantages: [
                        'No need for starting cell specification',
                        'Handles complex, branching trajectories',
                        'Scales to large datasets',
                        'Robust to noise'
                    ],
                    databases: [
                        { name: 'Monocle 3', url: 'https://cole-trapnell-lab.github.io/monocle3/', description: 'Official documentation' }
                    ]
                },
                {
                    id: 'scvelo',
                    name: 'scVelo',
                    description: 'RNA velocity analysis for cell fate prediction',
                    specifications: {
                        language: 'Python',
                        method: 'Spliced/unspliced RNA dynamics',
                        citation: 'Bergen et al., Nature Biotechnology 2020',
                        models: 'Deterministic, stochastic, dynamical'
                    },
                    features: [
                        'Velocity estimation from spliced/unspliced counts',
                        'Dynamic model fitting',
                        'Latent time recovery',
                        'Driver gene identification',
                        'Velocity confidence estimation',
                        'Stream plots and embeddings'
                    ],
                    advantages: [
                        'Reveals RNA dynamics',
                        'Predicts future cell states',
                        'Model-based approach',
                        'Identifies key regulatory genes'
                    ],
                    databases: [
                        { name: 'scVelo', url: 'https://scvelo.readthedocs.io/', description: 'Documentation and tutorials' },
                        { name: 'RNA Velocity', url: 'http://velocyto.org/', description: 'Original velocyto framework' }
                    ]
                },
                {
                    id: 'cellrank',
                    name: 'CellRank',
                    description: 'Probabilistic fate mapping combining RNA velocity and transcriptional similarity',
                    specifications: {
                        language: 'Python',
                        method: 'Markov chains on velocity field',
                        citation: 'Lange et al., Nature Methods 2022'
                    },
                    features: [
                        'Combines velocity with expression similarity',
                        'Identifies initial and terminal states',
                        'Computes fate probabilities',
                        'Gene expression trends along lineages',
                        'Circular trajectory support'
                    ],
                    databases: [
                        { name: 'CellRank', url: 'https://cellrank.readthedocs.io/', description: 'Documentation' }
                    ]
                },
                {
                    id: 'cellphonedb',
                    name: 'CellPhoneDB',
                    description: 'Cell-cell communication inference from single-cell transcriptomics',
                    specifications: {
                        language: 'Python',
                        database: '3,500+ curated interactions',
                        statistical: 'Permutation testing',
                        version: 'v3+ with spatial support'
                    },
                    features: [
                        'Ligand-receptor interaction prediction',
                        'Statistical significance testing',
                        'Multisubunit complex handling',
                        'Spatial interaction analysis',
                        'Deconvoluted data support'
                    ],
                    databases: [
                        { name: 'CellPhoneDB', url: 'https://www.cellphonedb.org/', description: 'Web interface and database' }
                    ]
                },
                {
                    id: 'nichenetr',
                    name: 'NicheNet',
                    description: 'Predict ligand-target gene regulatory potential',
                    specifications: {
                        language: 'R',
                        method: 'Prior knowledge integration',
                        database: 'Curated signaling networks'
                    },
                    features: [
                        'Ligand activity prediction',
                        'Target gene prediction',
                        'Ligand-receptor-target linking',
                        'Differential NicheNet analysis'
                    ],
                    databases: [
                        { name: 'NicheNet', url: 'https://github.com/saeyslab/nichenetr', description: 'Package and tutorials' }
                    ]
                },
                {
                    id: 'scenic',
                    name: 'SCENIC',
                    description: 'Single-cell regulatory network inference and clustering',
                    specifications: {
                        language: 'R/Python',
                        method: 'Coexpression + motif enrichment',
                        citation: 'Aibar et al., Nature Methods 2017'
                    },
                    versions: [
                        'SCENIC - Original implementation',
                        'pySCENIC - Python version',
                        'SCENIC+ - Chromatin accessibility integration'
                    ],
                    features: [
                        'Gene regulatory network inference',
                        'Regulon activity scoring',
                        'Cell state identification',
                        'Key transcription factor discovery'
                    ],
                    databases: [
                        { name: 'SCENIC', url: 'https://scenic.aertslab.org/', description: 'Documentation and resources' }
                    ]
                }
            ]
        },
        {
            id: 'visualization',
            name: 'Visualization Tools',
            icon: 'fas fa-eye',
            color: 'warning',
            items: [
                {
                    id: 'cellxgene',
                    name: 'cellxgene',
                    description: 'Interactive single-cell data explorer',
                    specifications: {
                        type: 'Web application',
                        format: 'h5ad (AnnData), h5 (loom)',
                        deployment: 'Local/Cloud/Hosted',
                        scalability: 'Millions of cells'
                    },
                    features: [
                        'Interactive UMAP/tSNE exploration',
                        'Gene expression overlay',
                        'Cell subset selection',
                        'Differential expression analysis',
                        'Annotation editing',
                        'Data hosting platform (CELLxGENE Discover)'
                    ],
                    deployment: [
                        'cellxgene launch data.h5ad - Local viewing',
                        'cellxgene prepare - Optimize for deployment',
                        'VIP extension - Additional features',
                        'Gateway - Multi-dataset hosting'
                    ],
                    databases: [
                        { name: 'cellxgene', url: 'https://cellxgene.cziscience.com/', description: 'Data portal and viewer' },
                        { name: 'cellxgene VIP', url: 'https://github.com/interactivereport/cellxgene_VIP', description: 'Extended features' }
                    ]
                },
                {
                    id: 'ucsc-cellbrowser',
                    name: 'UCSC Cell Browser',
                    description: 'Web-based visualization for single-cell datasets',
                    specifications: {
                        hosting: 'Self-hosted or UCSC hosted',
                        format: 'Multiple (h5ad, Seurat, text)',
                        public: 'cells.ucsc.edu'
                    },
                    features: [
                        'Dataset collection hosting',
                        'Gene expression visualization',
                        'Metadata filtering and coloring',
                        'Trajectory visualization',
                        'Dataset comparison',
                        'Publication-ready figures'
                    ],
                    databases: [
                        { name: 'UCSC Cell Browser', url: 'https://cells.ucsc.edu/', description: 'Public datasets portal' },
                        { name: 'Documentation', url: 'https://cellbrowser.readthedocs.io/', description: 'Setup and usage guide' }
                    ]
                },
                {
                    id: 'vitessce',
                    name: 'Vitessce',
                    description: 'Visual integration tool for exploration of spatial single-cell experiments',
                    specifications: {
                        type: 'JavaScript/React framework',
                        spatial: true,
                        multimodal: true,
                        embedding: 'Web-embeddable'
                    },
                    features: [
                        'Coordinated multiple views',
                        'Spatial data visualization',
                        'Gene expression heatmaps',
                        'Cell type hierarchies',
                        'Genomic profiles',
                        'Custom view configurations'
                    ],
                    databases: [
                        { name: 'Vitessce', url: 'http://vitessce.io/', description: 'Documentation and demos' }
                    ]
                },
                {
                    id: 'cirrocumulus',
                    name: 'Cirrocumulus',
                    description: 'Fast interactive visualization for millions of cells',
                    specifications: {
                        type: 'Web application',
                        scalability: '10+ million cells',
                        format: 'h5ad, zarr, Seurat'
                    },
                    features: [
                        'WebGL-accelerated rendering',
                        'Server-based computation',
                        'Spatial data support',
                        'Differential expression',
                        'Dataset sharing'
                    ],
                    databases: [
                        { name: 'Cirrocumulus', url: 'https://cirrocumulus.readthedocs.io/', description: 'Documentation' }
                    ]
                },
                {
                    id: 'scimilarity',
                    name: 'SCimilarity',
                    description: 'Web portal for queryable single-cell reference atlas',
                    specifications: {
                        type: 'Web application',
                        method: 'Deep learning embeddings',
                        reference: '22+ million cells'
                    },
                    features: [
                        'Cell type query and annotation',
                        'Similar cell finding',
                        'Reference atlas exploration',
                        'Batch effect robust'
                    ],
                    databases: [
                        { name: 'SCimilarity', url: 'https://www.scimilarity.org/', description: 'Web portal' }
                    ]
                }
            ]
        },
        {
            id: 'integration-tools',
            name: 'Data Integration',
            icon: 'fas fa-link',
            color: 'info',
            items: [
                {
                    id: 'harmony',
                    name: 'Harmony',
                    description: 'Fast batch correction for single-cell data',
                    specifications: {
                        language: 'R/Python',
                        method: 'Iterative clustering correction',
                        speed: 'Linear scaling',
                        citation: 'Korsunsky et al., Nature Methods 2019'
                    },
                    features: [
                        'Fast batch correction',
                        'Preserves biological variation',
                        'Works in reduced dimensions',
                        'Multiple batch variables'
                    ],
                    databases: [
                        { name: 'Harmony', url: 'https://github.com/immunogenomics/harmony', description: 'Package and tutorials' }
                    ]
                },
                {
                    id: 'scanorama',
                    name: 'Scanorama',
                    description: 'Panoramic stitching of heterogeneous single-cell data',
                    specifications: {
                        language: 'Python',
                        method: 'Mutual nearest neighbors',
                        scalability: 'Millions of cells'
                    },
                    features: [
                        'Batch correction',
                        'Cross-dataset projection',
                        'Feature matching',
                        'Handles heterogeneous datasets'
                    ],
                    databases: [
                        { name: 'Scanorama', url: 'https://github.com/brianhie/scanorama', description: 'Implementation' }
                    ]
                },
                {
                    id: 'liger',
                    name: 'LIGER',
                    description: 'Integrative non-negative matrix factorization for single-cell data',
                    specifications: {
                        language: 'R',
                        method: 'iNMF',
                        multimodal: true,
                        citation: 'Welch et al., Cell 2019'
                    },
                    features: [
                        'Cross-species integration',
                        'Multimodal integration',
                        'Shared and specific factors',
                        'Developmental trajectory alignment'
                    ],
                    databases: [
                        { name: 'LIGER', url: 'https://github.com/welch-lab/liger', description: 'Package documentation' }
                    ]
                },
                {
                    id: 'scib',
                    name: 'scib',
                    description: 'Benchmarking atlas-level data integration',
                    specifications: {
                        language: 'Python',
                        metrics: '14 integration metrics',
                        methods: '16+ integration methods'
                    },
                    features: [
                        'Comprehensive benchmarking',
                        'Biological conservation metrics',
                        'Batch correction metrics',
                        'Method comparison'
                    ],
                    databases: [
                        { name: 'scib', url: 'https://scib.readthedocs.io/', description: 'Benchmarking framework' }
                    ]
                }
            ]
        },
        {
            id: 'annotation-tools',
            name: 'Cell Type Annotation',
            icon: 'fas fa-tags',
            color: 'danger',
            items: [
                {
                    id: 'singleR',
                    name: 'SingleR',
                    description: 'Reference-based automatic cell type annotation',
                    specifications: {
                        language: 'R',
                        method: 'Correlation-based',
                        references: 'Multiple built-in'
                    },
                    features: [
                        'Automatic annotation',
                        'Multiple reference support',
                        'Fine-tuning option',
                        'Confidence scores'
                    ],
                    databases: [
                        { name: 'SingleR', url: 'https://bioconductor.org/packages/SingleR/', description: 'Bioconductor package' }
                    ]
                },
                {
                    id: 'celltypist',
                    name: 'CellTypist',
                    description: 'Automated cell type annotation with pre-trained models',
                    specifications: {
                        language: 'Python',
                        models: '50+ pre-trained',
                        method: 'Logistic regression'
                    },
                    features: [
                        'Multiple tissue models',
                        'Cross-species annotation',
                        'Model training capability',
                        'Majority voting'
                    ],
                    databases: [
                        { name: 'CellTypist', url: 'https://www.celltypist.org/', description: 'Models and documentation' }
                    ]
                },
                {
                    id: 'sctype',
                    name: 'scType',
                    description: 'Marker gene-based cell type annotation',
                    specifications: {
                        language: 'R',
                        database: 'Cell marker database',
                        method: 'Gene set scoring'
                    },
                    features: [
                        'Marker gene database',
                        'Custom marker support',
                        'Hierarchical annotation',
                        'Visualization tools'
                    ],
                    databases: [
                        { name: 'scType', url: 'https://github.com/IanevskiAleksandr/sc-type', description: 'Package and database' }
                    ]
                },
                {
                    id: 'azimuth',
                    name: 'Azimuth',
                    description: 'Web app for reference-based single-cell analysis',
                    specifications: {
                        type: 'Web application',
                        references: 'Human and mouse atlases',
                        method: 'Seurat mapping'
                    },
                    features: [
                        'Web-based annotation',
                        'Multiple reference atlases',
                        'UMAP projection',
                        'Downloadable results'
                    ],
                    databases: [
                        { name: 'Azimuth', url: 'https://azimuth.hubmapconsortium.org/', description: 'Web application' }
                    ]
                }
            ]
        }
    ],
    
    resources: {
        databases: [
            { name: 'scRNA-tools', url: 'https://www.scrna-tools.org/', type: 'Comprehensive database of >1000 scRNA-seq analysis tools' },
            { name: 'Awesome Single Cell', url: 'https://github.com/seandavi/awesome-single-cell', type: 'Curated list of single-cell analysis methods' },
            { name: 'Single Cell Portal', url: 'https://singlecell.broadinstitute.org/', type: 'Data repository with analysis tools' },
            { name: 'CELLxGENE Discover', url: 'https://cellxgene.cziscience.com/', type: 'Standardized single-cell data repository' },
            { name: 'ASAP', url: 'https://asap.epfl.ch/', type: 'Automated Single-cell Analysis Portal' },
            { name: 'scRNAseq Benchmark', url: 'https://github.com/tabdelaal/scRNAseq_Benchmark', type: 'Systematic evaluation of methods' }
        ],
        tools: [
            { name: 'Galaxy Single Cell', url: 'https://singlecell.usegalaxy.eu/', purpose: 'Web-based analysis platform without coding' },
            { name: 'DISCO', url: 'https://www.immunesinglecell.org/', purpose: 'Database of Immune Single Cell Expression' },
            { name: 'SCPortalen', url: 'https://single-cell.clst.riken.jp/', purpose: 'RIKEN single-cell database' },
            { name: 'PanglaoDB', url: 'https://panglaodb.se/', purpose: 'Single-cell RNA sequencing database' },
            { name: 'CellMarker 2.0', url: 'http://bio-bigdata.hrbmu.edu.cn/CellMarker/', purpose: 'Cell marker database for annotation' }
        ],
        tutorials: [
            { title: 'OSCA Book', url: 'https://bioconductor.org/books/release/OSCA/', description: 'Orchestrating Single-Cell Analysis with Bioconductor' },
            { title: 'Best Practices', url: 'https://www.sc-best-practices.org/', description: 'Current best practices in single-cell analysis' },
            { title: 'Seurat Vignettes', url: 'https://satijalab.org/seurat/articles/get_started.html', description: 'Comprehensive Seurat tutorials' },
            { title: 'Scanpy Tutorials', url: 'https://scanpy-tutorials.readthedocs.io/', description: 'Collection of Scanpy analysis workflows' },
            { title: 'HCA Methods', url: 'https://www.humancellatlas.org/analyze/', description: 'Human Cell Atlas analysis methods' }
        ]
    }
};