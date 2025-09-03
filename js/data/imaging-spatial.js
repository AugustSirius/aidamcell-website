export const imagingSpatial = {
    id: 'imaging-spatial',
    name: 'Imaging & Spatial Technologies',
    description: 'Advanced microscopy, spatial transcriptomics, and imaging technologies for single-cell and tissue analysis',
    
    categories: [
        {
            id: 'spatial-sequencing',
            name: 'Spatial Sequencing Technologies',
            icon: 'fas fa-map-marked-alt',
            color: 'primary',
            items: [
                {
                    id: 'visium',
                    name: '10x Genomics Visium',
                    description: 'Spatial gene expression profiling with morphological context preservation',
                    specifications: {
                        resolution: '55 µm spot diameter',
                        spotsPerSlide: '~5,000 spots',
                        genes: 'Whole transcriptome (~18,000 genes detected)',
                        tissueArea: '6.5 x 6.5 mm capture area',
                        compatibility: 'Fresh frozen, FFPE, Fixed frozen'
                    },
                    workflow: [
                        'Tissue sectioning (10 µm thickness)',
                        'H&E or immunofluorescence staining',
                        'Brightfield or fluorescence imaging',
                        'Tissue permeabilization (12-24 minutes)',
                        'Reverse transcription on-slide',
                        'Second strand synthesis and denaturation',
                        'cDNA amplification and library construction',
                        'Sequencing on Illumina platforms'
                    ],
                    advantages: [
                        'Unbiased whole transcriptome coverage',
                        'Direct integration with histology',
                        'Established commercial platform',
                        'Compatible with various tissue types',
                        'Strong computational support'
                    ],
                    limitations: [
                        'Not single-cell resolution (2-10 cells per spot)',
                        'Limited to thin tissue sections',
                        'High cost per sample (~$1,500-3,000)',
                        'Requires optimization for each tissue type'
                    ],
                    databases: [
                        { name: '10x Genomics Datasets', url: 'https://www.10xgenomics.com/resources/datasets', description: 'Official Visium example datasets' },
                        { name: 'SpatialDB', url: 'https://www.spatialomics.org/SpatialDB/', description: 'Spatial transcriptomics database' },
                        { name: 'STOmicsDB', url: 'https://db.cngb.org/stomics/', description: 'Spatiotemporal omics database' }
                    ]
                },
                {
                    id: 'visium-hd',
                    name: 'Visium HD',
                    description: 'High-resolution spatial gene expression from 10x Genomics',
                    specifications: {
                        resolution: '2 µm x 2 µm bins',
                        continuous: 'Continuous squares, no gaps',
                        genes: 'Whole transcriptome',
                        tissueArea: '6.5 x 6.5 mm',
                        flexibility: 'Multiple bin sizes (2, 8, 16 µm)'
                    },
                    advantages: [
                        'Near single-cell resolution',
                        'Flexible binning post-capture',
                        'No spot boundaries',
                        'Better cellular resolution'
                    ],
                    limitations: [
                        'Higher sequencing requirements',
                        'More complex analysis',
                        'Limited availability (2024 launch)'
                    ]
                },
                {
                    id: 'slide-seq',
                    name: 'Slide-seq/Slide-seqV2',
                    description: 'High-resolution spatial transcriptomics using DNA-barcoded beads',
                    specifications: {
                        resolution: '10 µm bead diameter',
                        coverage: 'Whole transcriptome',
                        beadDensity: '~10,000 beads/mm²',
                        tissueSize: 'Up to 10 mm diameter pucks'
                    },
                    workflow: [
                        'Array generation with barcoded beads',
                        'Tissue section placement on array',
                        'In situ reverse transcription',
                        'Tissue digestion and library prep',
                        'Sequencing and spatial reconstruction'
                    ],
                    advantages: [
                        'Higher resolution than Visium',
                        'Unbiased transcriptome capture',
                        'Large tissue coverage',
                        'Open-source protocol'
                    ],
                    limitations: [
                        'Complex array preparation',
                        'Lower capture efficiency than Visium',
                        'Requires specialized equipment',
                        'No commercial kit available'
                    ],
                    databases: [
                        { name: 'Broad Institute Single Cell Portal', url: 'https://singlecell.broadinstitute.org/single_cell', description: 'Includes Slide-seq datasets' }
                    ]
                },
                {
                    id: 'stereo-seq',
                    name: 'Stereo-seq (STOmics)',
                    description: 'Ultra-large field nanoscale spatial omics platform from BGI',
                    specifications: {
                        resolution: '220-715 nm bin size',
                        chipSize: '13.2 x 13.2 cm',
                        spotsPerChip: 'Up to 20 billion DNB spots',
                        genes: 'Whole transcriptome'
                    },
                    advantages: [
                        'Largest capture area available',
                        'Subcellular resolution possible',
                        'Whole organ/embryo mapping',
                        'Cost-effective for large tissues'
                    ],
                    limitations: [
                        'Limited availability outside China',
                        'Requires BGI sequencing platform',
                        'Large data management requirements'
                    ],
                    applications: [
                        'Mouse embryo spatiotemporal atlas',
                        'Whole organ mapping',
                        'Large tissue sections',
                        'Developmental biology'
                    ],
                    databases: [
                        { name: 'MOSTA', url: 'https://db.cngb.org/stomics/mosta/', description: 'Mouse Organogenesis Spatiotemporal Transcriptomic Atlas' }
                    ]
                },
                {
                    id: 'dbit-seq',
                    name: 'DBiT-seq',
                    description: 'Deterministic Barcoding in Tissue for spatial omics',
                    specifications: {
                        resolution: '10-50 µm',
                        modalities: 'RNA, protein, or both',
                        barcoding: 'Microfluidic delivery',
                        coverage: 'Transcriptome or targeted panels'
                    },
                    advantages: [
                        'Multi-omics capability',
                        'No specialized slides needed',
                        'Works with standard equipment',
                        'Flexible resolution'
                    ],
                    workflow: [
                        'Tissue fixation and permeabilization',
                        'First set of DNA barcodes (rows)',
                        'Second set of DNA barcodes (columns)',
                        'In situ reverse transcription',
                        'Library preparation and sequencing'
                    ]
                }
            ]
        },
        {
            id: 'imaging-based-spatial',
            name: 'Imaging-Based Spatial Methods',
            icon: 'fas fa-images',
            color: 'secondary',
            items: [
                {
                    id: 'merfish',
                    name: 'MERFISH',
                    description: 'Multiplexed error-robust fluorescence in situ hybridization',
                    specifications: {
                        genes: '100-10,000 genes',
                        resolution: 'Subcellular (~100 nm)',
                        accuracy: '>95% with error correction',
                        zStack: '10-20 µm tissue depth',
                        imagingRounds: '8-20 rounds'
                    },
                    principle: 'Combinatorial labeling with binary barcoding and error correction codes',
                    workflow: [
                        'Probe library design with error-correcting codes',
                        'Sample fixation and permeabilization',
                        'Primary probe hybridization',
                        'Sequential rounds of secondary probe imaging',
                        'Image registration and alignment',
                        'Spot detection and barcode decoding',
                        'Cell segmentation and quantification'
                    ],
                    advantages: [
                        'True single-cell resolution',
                        'Subcellular RNA localization',
                        '3D tissue reconstruction capability',
                        'High detection accuracy',
                        'Compatible with protein co-detection'
                    ],
                    limitations: [
                        'Limited to selected gene panels',
                        'Long imaging time (24-48 hours)',
                        'Expensive probe synthesis',
                        'Requires specialized microscope'
                    ],
                    commercialPlatforms: [
                        'Vizgen MERSCOPE',
                        'Rebus Esper'
                    ],
                    databases: [
                        { name: 'Vizgen Data Release', url: 'https://info.vizgen.com/mouse-brain-data-release', description: 'MERFISH mouse brain atlas' },
                        { name: 'Allen Brain Map', url: 'https://portal.brain-map.org/', description: 'Includes MERFISH datasets' }
                    ]
                },
                {
                    id: 'xenium',
                    name: '10x Genomics Xenium',
                    description: 'In situ gene expression analysis at subcellular resolution',
                    specifications: {
                        genes: 'Up to 5,000 custom genes',
                        resolution: 'Subcellular (~200 nm)',
                        tissueSize: 'Up to 2.8 x 1.6 cm',
                        workflow: 'Automated, 24-48 hours',
                        zDepth: '10-20 µm sections'
                    },
                    features: [
                        'Automated end-to-end workflow',
                        'Pre-designed and custom panels',
                        'Multimodal readouts (RNA + protein)',
                        'H&E staining compatible',
                        'Advanced cell segmentation'
                    ],
                    panels: {
                        'Human': 'Breast, lung, brain, kidney, skin',
                        'Mouse': 'Brain, tissue multi-panel',
                        'Custom': 'User-defined up to 480 genes'
                    },
                    advantages: [
                        'User-friendly instrument',
                        'High sensitivity and specificity',
                        'Robust cell segmentation',
                        'Commercial support'
                    ],
                    databases: [
                        { name: '10x Xenium Explorer', url: 'https://www.10xgenomics.com/products/xenium-explorer', description: 'Visualization software and datasets' }
                    ]
                },
                {
                    id: 'seqfish',
                    name: 'seqFISH+',
                    description: 'Sequential fluorescence in situ hybridization with pseudo-color encoding',
                    specifications: {
                        genes: 'Up to 10,000 genes',
                        resolution: 'Subcellular',
                        rounds: '80 rounds of imaging',
                        accuracy: '>80% detection efficiency'
                    },
                    methodology: 'Pseudo-color barcoding with sequential rounds',
                    advantages: [
                        'Highest gene count for imaging methods',
                        'Single molecule detection',
                        'Compatible with immunofluorescence',
                        'Preserves tissue architecture'
                    ],
                    limitations: [
                        'Extremely long protocol (3-4 days)',
                        'Complex probe design',
                        'Photobleaching concerns',
                        'Limited to research settings'
                    ],
                    keyPublications: [
                        'Eng et al., Nature 2019 - 10,000 genes',
                        'Shah et al., Cell 2018 - intron seqFISH'
                    ]
                },
                {
                    id: 'cosmx',
                    name: 'NanoString CosMx SMI',
                    description: 'Spatial molecular imager for RNA and protein',
                    specifications: {
                        rna: '1,000+ transcripts',
                        protein: '64+ proteins',
                        resolution: 'Subcellular',
                        fov: 'Multiple fields of view',
                        tissueCompatibility: 'FFPE, fresh frozen'
                    },
                    technology: 'Cyclic in situ hybridization with reporter probes',
                    workflow: [
                        'Tissue preparation and imaging',
                        'Target probe hybridization',
                        'Cyclic reporter hybridization',
                        'UV cleavage between cycles',
                        'Image processing and analysis'
                    ],
                    advantages: [
                        'RNA and protein co-detection',
                        'FFPE compatibility',
                        'High-plex capability',
                        'Automated instrument'
                    ],
                    applications: [
                        'Tumor microenvironment',
                        'Immunology studies',
                        'Neuroscience',
                        'Cell atlas projects'
                    ],
                    databases: [
                        { name: 'NanoString Spatial Organ Atlas', url: 'https://nanostring.com/products/cosmx-spatial-molecular-imager/spatial-organ-atlas/', description: 'Reference datasets' }
                    ]
                },
                {
                    id: 'iss',
                    name: 'In Situ Sequencing',
                    description: 'Direct sequencing of RNA in tissue sections',
                    specifications: {
                        targets: '50-500 genes',
                        resolution: 'Subcellular',
                        accuracy: '99% base calling',
                        readLength: '4-5 bases typical'
                    },
                    methods: [
                        'Padlock probes and RCA',
                        'Sequencing by ligation',
                        'HybISS (gap-filling approach)',
                        'FISSEQ (fluorescent sequencing)'
                    ],
                    advantages: [
                        'Direct sequence readout',
                        'High specificity',
                        'Mutation detection capability',
                        'Isoform detection possible'
                    ],
                    limitations: [
                        'Limited read length',
                        'Lower throughput',
                        'Complex chemistry'
                    ]
                }
            ]
        },
        {
            id: 'mass-spec-imaging',
            name: 'Mass Spectrometry Imaging',
            icon: 'fas fa-atom',
            color: 'success',
            items: [
                {
                    id: 'imaging-mass-cytometry',
                    name: 'Imaging Mass Cytometry (IMC)',
                    description: 'Highly multiplexed protein imaging using metal-tagged antibodies',
                    specifications: {
                        markers: '40+ simultaneous proteins',
                        resolution: '1 µm',
                        tissueArea: 'Up to 1 mm²',
                        depth: 'Surface ablation (~200 nm)',
                        quantitative: true
                    },
                    principle: 'Laser ablation of tissue with time-of-flight mass spectrometry detection',
                    workflow: [
                        'Metal-conjugated antibody staining',
                        'Laser ablation pixel by pixel',
                        'Mass cytometry detection',
                        'Image reconstruction',
                        'Cell segmentation and analysis'
                    ],
                    advantages: [
                        'No spectral overlap',
                        'Highly quantitative',
                        'Works with FFPE',
                        'Stable metal reporters'
                    ],
                    limitations: [
                        'Destructive imaging',
                        'Slow acquisition (1 hr/mm²)',
                        'Expensive antibodies',
                        'Limited to surface proteins'
                    ],
                    vendor: 'Standard BioTools (formerly Fluidigm)',
                    databases: [
                        { name: 'MIBI-TOF Data', url: 'https://www.angelolab.com/mibi-data', description: 'Example IMC/MIBI datasets' }
                    ]
                },
                {
                    id: 'mibi-tof',
                    name: 'MIBI-TOF',
                    description: 'Multiplexed ion beam imaging by time of flight',
                    specifications: {
                        markers: '40+ proteins',
                        resolution: '260 nm',
                        speed: 'Faster than IMC',
                        sensitivity: 'Single molecule'
                    },
                    technology: 'Secondary ion mass spectrometry (SIMS)',
                    advantages: [
                        'Superior spatial resolution',
                        'Subcellular features visible',
                        'Highly sensitive',
                        'Quantitative measurements'
                    ],
                    applications: [
                        'Tumor immune microenvironment',
                        'Subcellular protein localization',
                        'Tissue architecture studies'
                    ],
                    vendor: 'Ionpath',
                    databases: [
                        { name: 'MIBI Tumor Atlas', url: 'https://www.tissue-atlas.org/', description: 'Tumor microenvironment atlas' }
                    ]
                },
                {
                    id: 'maldi-imaging',
                    name: 'MALDI Imaging',
                    description: 'Matrix-assisted laser desorption/ionization for spatial metabolomics',
                    specifications: {
                        resolution: '5-50 µm',
                        molecules: 'Metabolites, lipids, peptides, drugs',
                        massRange: '50-30,000 Da',
                        tissueSize: 'Whole tissue sections'
                    },
                    advantages: [
                        'Label-free detection',
                        'Broad molecular coverage',
                        'Drug distribution studies',
                        'Metabolite mapping'
                    ],
                    limitations: [
                        'Matrix application variability',
                        'Ion suppression effects',
                        'Limited protein detection',
                        'Requires careful sample prep'
                    ],
                    applications: [
                        'Drug ADME studies',
                        'Biomarker discovery',
                        'Lipidomics',
                        'Cancer metabolism'
                    ]
                },
                {
                    id: 'desi-imaging',
                    name: 'DESI-MS Imaging',
                    description: 'Desorption electrospray ionization mass spectrometry imaging',
                    specifications: {
                        resolution: '35-200 µm',
                        speed: 'Fast (5-10 pixels/second)',
                        ambient: 'Minimal sample preparation',
                        molecules: 'Lipids, metabolites, drugs'
                    },
                    advantages: [
                        'Ambient conditions',
                        'Minimal sample prep',
                        'Fast acquisition',
                        'Soft ionization'
                    ],
                    applications: [
                        'Surgical margin assessment',
                        'Metabolomics',
                        'Drug imaging',
                        'Biomarker discovery'
                    ]
                }
            ]
        },
        {
            id: 'multiplexed-if',
            name: 'Multiplexed Immunofluorescence',
            icon: 'fas fa-palette',
            color: 'warning',
            items: [
                {
                    id: 'codex',
                    name: 'CODEX/PhenoCycler',
                    description: 'CO-Detection by indEXing for highly multiplexed tissue imaging',
                    specifications: {
                        markers: '40-60+ proteins',
                        resolution: 'Single-cell in tissue',
                        cycles: '15-20 cycles',
                        tissueType: 'FFPE, fresh frozen',
                        automation: 'Fluidics system available'
                    },
                    principle: 'DNA-barcoded antibodies with cyclic fluorescent reporter addition',
                    workflow: [
                        'Antibody-oligonucleotide conjugation',
                        'Simultaneous staining of all targets',
                        'Cyclic imaging with fluorescent oligos',
                        'Gentle stripping between cycles',
                        'Image registration and processing'
                    ],
                    advantages: [
                        'Works with standard fluorescence microscopes',
                        'All antibodies applied once',
                        'Gentle on tissue morphology',
                        '3D imaging capability'
                    ],
                    limitations: [
                        'Long acquisition time',
                        'Antibody validation needed',
                        'Autofluorescence management'
                    ],
                    vendor: 'Akoya Biosciences (PhenoCycler)',
                    databases: [
                        { name: 'HuBMAP', url: 'https://portal.hubmapconsortium.org/', description: 'Human BioMolecular Atlas includes CODEX data' }
                    ]
                },
                {
                    id: 'cycif',
                    name: 't-CyCIF',
                    description: 'Tissue-based cyclic immunofluorescence',
                    specifications: {
                        markers: '60+ proteins',
                        cycles: '10-25 rounds',
                        resolution: 'Subcellular',
                        tissueCompatibility: 'FFPE optimal'
                    },
                    methodology: 'Cyclic antibody staining, imaging, and removal',
                    advantages: [
                        'Uses conventional antibodies',
                        'No specialized equipment',
                        'Open-source protocols',
                        'High-quality images'
                    ],
                    workflow: [
                        'FFPE section preparation',
                        'Antibody staining (3-4 markers)',
                        'Imaging',
                        'Antibody removal (bleaching)',
                        'Repeat cycles'
                    ],
                    databases: [
                        { name: 'CyCIF.org', url: 'https://www.cycif.org/', description: 'Methods and example data' }
                    ]
                },
                {
                    id: 'ibex',
                    name: 'IBEX',
                    description: 'Iterative bleaching extends multiplexity',
                    specifications: {
                        markers: '65+ proteins demonstrated',
                        bleaching: 'Lithium borohydride',
                        compatibility: 'Over 270 validated antibodies',
                        openSource: true
                    },
                    advantages: [
                        'Simple chemistry',
                        'Works with diverse tissues',
                        'Community antibody database',
                        'Cost-effective'
                    ],
                    resources: [
                        { name: 'IBEX Knowledge Base', url: 'https://github.com/IBEXImagingCommunity', description: 'Protocols and antibody database' }
                    ]
                },
                {
                    id: '4i',
                    name: '4i (Iterative Indirect Immunofluorescence)',
                    description: 'Multiplexed protein imaging through iterative cycles',
                    specifications: {
                        markers: '40-100 proteins',
                        resolution: 'Single-cell',
                        elution: 'Chemical elution',
                        automation: 'Robot-assisted'
                    },
                    advantages: [
                        'Gentle elution conditions',
                        'Maintains morphology',
                        'Automated workflows available',
                        'Whole slide imaging'
                    ],
                    applications: [
                        'Drug response profiling',
                        'Tumor heterogeneity',
                        'Signaling pathway analysis'
                    ]
                }
            ]
        },
        {
            id: 'advanced-microscopy',
            name: 'Advanced Microscopy Techniques',
            icon: 'fas fa-microscope',
            color: 'info',
            items: [
                {
                    id: 'light-sheet',
                    name: 'Light Sheet Fluorescence Microscopy',
                    description: 'Fast 3D imaging with minimal photodamage for live samples',
                    specifications: {
                        speed: '100-1000 fps',
                        resolution: 'XY: 250 nm, Z: 1-2 µm',
                        sampleSize: 'Cells to whole organs',
                        duration: 'Hours to days of imaging'
                    },
                    variants: [
                        'SPIM (Selective Plane Illumination)',
                        'Lattice light sheet',
                        'Swept confocally aligned planar excitation (SCAPE)',
                        'Oblique plane microscopy (OPM)'
                    ],
                    advantages: [
                        'Minimal phototoxicity',
                        'Fast volumetric imaging',
                        'Long-term live imaging',
                        'Large sample capability'
                    ],
                    applications: [
                        'Embryo development',
                        'Organoid imaging',
                        'Cleared tissue imaging',
                        'Cell migration tracking'
                    ],
                    systems: [
                        'Zeiss Lightsheet 7',
                        'Bruker Luxendo',
                        'Miltenyi UltraMicroscope'
                    ]
                },
                {
                    id: 'super-resolution',
                    name: 'Super-Resolution Microscopy',
                    description: 'Breaking the diffraction limit for nanoscale imaging',
                    techniques: {
                        'STORM/PALM': '20-30 nm resolution, single molecule localization',
                        'STED': '30-80 nm, stimulated emission depletion',
                        'SIM': '100-130 nm, structured illumination',
                        'MINFLUX': '<5 nm, minimal photon fluxes',
                        'Expansion microscopy': 'Physical sample expansion 4-20x'
                    },
                    applications: [
                        'Protein complex architecture',
                        'Synaptic structure',
                        'Chromatin organization',
                        'Membrane domains'
                    ],
                    advantages: [
                        'Molecular-scale resolution',
                        'Multi-color capability',
                        'Live cell compatible (some methods)',
                        'Quantitative measurements'
                    ],
                    limitations: [
                        'Limited imaging depth',
                        'Specialized equipment',
                        'Complex sample preparation',
                        'Photobleaching concerns'
                    ]
                },
                {
                    id: 'two-photon',
                    name: 'Two-Photon Microscopy',
                    description: 'Deep tissue imaging in living samples using nonlinear excitation',
                    specifications: {
                        penetration: 'Up to 1 mm depth',
                        resolution: 'XY: 0.5 µm, Z: 2 µm',
                        wavelength: '700-1300 nm excitation',
                        applications: 'Brain, intravital imaging'
                    },
                    advantages: [
                        'Deep tissue penetration',
                        'Reduced photodamage',
                        'Intrinsic optical sectioning',
                        'Less scattering'
                    ],
                    variants: [
                        'Three-photon microscopy (deeper)',
                        'Adaptive optics two-photon',
                        'Resonant scanning',
                        'Random access multiphoton'
                    ],
                    databases: [
                        { name: 'Allen Brain Observatory', url: 'https://observatory.brain-map.org/', description: 'Two-photon calcium imaging data' }
                    ]
                },
                {
                    id: 'cryo-em',
                    name: 'Cryo-Electron Microscopy',
                    description: 'Native structure imaging at near-atomic resolution',
                    specifications: {
                        resolution: '1-4 Å for single particles',
                        samplePrep: 'Vitrification',
                        temperature: '-196°C',
                        techniques: 'Single particle, tomography, MicroED'
                    },
                    applications: [
                        'Protein structure determination',
                        'Cellular ultrastructure',
                        'Virus structure',
                        'Membrane proteins'
                    ],
                    advantages: [
                        'Near-native state',
                        'No crystallization needed',
                        'Large complexes possible',
                        'Time-resolved studies'
                    ],
                    databases: [
                        { name: 'EMDB', url: 'https://www.ebi.ac.uk/emdb/', description: 'Electron Microscopy Data Bank' },
                        { name: 'EMPIAR', url: 'https://www.ebi.ac.uk/empiar/', description: 'EM Public Image Archive' }
                    ]
                },
                {
                    id: 'label-free-imaging',
                    name: 'Label-Free Imaging',
                    description: 'Visualization without fluorescent labels',
                    techniques: {
                        'Phase contrast': 'Refractive index differences',
                        'DIC': 'Differential interference contrast',
                        'Quantitative phase': 'Dry mass measurements',
                        'Raman microscopy': 'Molecular fingerprinting',
                        'SHG/THG': 'Second/third harmonic generation',
                        'OCT': 'Optical coherence tomography'
                    },
                    advantages: [
                        'No photobleaching',
                        'Minimal sample preparation',
                        'Live cell compatible',
                        'Quantitative measurements'
                    ],
                    applications: [
                        'Cell mass measurements',
                        'Organelle dynamics',
                        'Drug response',
                        'Tissue structure'
                    ]
                }
            ]
        }
    ],
    
    resources: {
        databases: [
            { name: 'Image Data Resource', url: 'https://idr.openmicroscopy.org/', type: 'Public microscopy datasets' },
            { name: 'Cell Image Library', url: 'http://www.cellimagelibrary.org/', type: 'Cellular imaging database' },
            { name: 'BioImage Archive', url: 'https://www.ebi.ac.uk/bioimage-archive/', type: 'Biological imaging data' },
            { name: 'SSBD', url: 'http://ssbd.qbic.riken.jp/', type: 'Systems Science of Biological Dynamics database' },
            { name: 'SpatialData', url: 'https://spatialdata.scverse.org/', type: 'Spatial omics data standard' }
        ],
        tools: [
            { name: 'QuPath', url: 'https://qupath.github.io/', purpose: 'Bioimage analysis for digital pathology' },
            { name: 'napari', url: 'https://napari.org/', purpose: 'Multi-dimensional image viewer' },
            { name: 'ImageJ/Fiji', url: 'https://fiji.sc/', purpose: 'General image analysis' },
            { name: 'CellProfiler', url: 'https://cellprofiler.org/', purpose: 'Cell image analysis' },
            { name: 'Squidpy', url: 'https://squidpy.readthedocs.io/', purpose: 'Spatial single-cell analysis' },
            { name: 'Giotto', url: 'https://giottosuite.readthedocs.io/', purpose: 'Spatial data analysis framework' },
            { name: 'MCMICRO', url: 'https://mcmicro.org/', purpose: 'Multiple-choice microscopy pipeline' },
            { name: 'Starfish', url: 'https://spacetx-starfish.readthedocs.io/', purpose: 'Image-based transcriptomics' }
        ],
        reviews: [
            { title: 'Museum of spatial transcriptomics', year: 2022, journal: 'Nature Methods', doi: '10.1038/s41592-022-01409-2' },
            { title: 'Spatial omics technologies', year: 2024, journal: 'Nature Reviews Methods Primers', doi: '10.1038/s43586-023-00276-1' },
            { title: 'The emerging landscape of spatial profiling', year: 2022, journal: 'Nature Reviews Genetics', doi: '10.1038/s41576-022-00515-3' }
        ]
    }
};