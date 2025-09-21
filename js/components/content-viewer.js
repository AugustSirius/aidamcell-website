export class ContentViewer {
    constructor(app) {
        this.app = app;
        this.container = document.getElementById('content-viewer');
    }
    
    render(sectionData, itemId = null) {
        if (!sectionData) {
            this.renderEmpty();
            return;
        }
        
        if (itemId) {
            this.renderDetail(sectionData, itemId);
        } else if (sectionData.id === 'overview') {
            this.renderOverview(sectionData);
        } else {
            this.renderSection(sectionData);
        }
    }
    
    renderOverview(data) {
        this.container.innerHTML = `
            <div class="overview">
                <h1>${data.name}</h1>
                <p class="description">${data.description}</p>
                
                <div class="stats-bar">
                    <div class="stat">
                        <div class="stat-value">${data.stats.totalItems}</div>
                        <div class="stat-label">Total Items</div>
                    </div>
                    <div class="stat">
                        <div class="stat-value">${data.stats.databases}</div>
                        <div class="stat-label">Databases</div>
                    </div>
                    <div class="stat">
                        <div class="stat-value">${data.stats.tools}</div>
                        <div class="stat-label">Tools</div>
                    </div>
                </div>
                
                <div class="section-grid">
                    ${data.highlights.map(section => this.createHighlightCard(section)).join('')}
                </div>
            </div>
        `;
    }
    
    renderSection(sectionData) {
        this.container.innerHTML = `
            <div class="section">
                <h1>${sectionData.name}</h1>
                <p class="description">${sectionData.description}</p>
                
                ${sectionData.categories.map(category => `
                    <div class="category">
                        <h2 class="category-title">
                            <i class="${category.icon}"></i>
                            ${category.name}
                        </h2>
                        <div class="section-grid">
                            ${category.items.map(item => this.createItemCard(item, sectionData.id)).join('')}
                        </div>
                    </div>
                `).join('')}
                
                ${this.renderResources(sectionData.resources)}
            </div>
        `;
        
        // Add click handlers
        this.attachCardHandlers();
    }
    
    renderDetail(sectionData, itemId) {
        const item = this.app.knowledgeBase.getItem(sectionData.id, itemId);
        if (!item) {
            this.renderEmpty();
            return;
        }
        
        this.container.innerHTML = `
            <div class="detail-container">
                <aside class="sidebar">
                    ${this.renderSidebar(sectionData, itemId)}
                </aside>
                
                <div class="main-content">
                    <h1>${item.name}</h1>
                    <p class="description">${item.description}</p>
                    
                    ${this.renderItemDetails(item)}
                </div>
            </div>
        `;
    }
    
    createHighlightCard(section) {
        return `
            <div class="card highlight-card" onclick="app.navigateToSection('${section.title.toLowerCase().replace(/ & /g, '-').replace(/ /g, '-')}')">
                <div class="card-header">
                    <div class="card-icon ${section.color}">
                        <i class="${section.icon}"></i>
                    </div>
                    <div>
                        <div class="card-title">${section.title}</div>
                        <div class="card-count">${section.count} categories</div>
                    </div>
                </div>
            </div>
        `;
    }
    
    createItemCard(item, sectionId) {
        // Pick the most suitable field for tags
        let rawTags = item.tags ?? item.components ?? item.platforms ?? item.methods ?? [];
    
        // Normalize to an array of strings
        if (Array.isArray(rawTags)) {
            // ok
        } else if (typeof rawTags === 'object' && rawTags !== null) {
            rawTags = Object.keys(rawTags);
        } else if (typeof rawTags === 'string') {
            rawTags = [rawTags];
        } else {
            rawTags = [];
        }
    
        return `
            <div class="card item-card"
                 data-section="${sectionId}"
                 data-item="${item.id}">
                <div class="card-header">
                    <div class="card-title">${item.name}</div>
                </div>
                <div class="card-description">${item.description}</div>
                <div class="card-tags">
                    ${rawTags.slice(0, 3)
                             .map(tag => `<span class="tag">${tag}</span>`)
                             .join('')}
                </div>
            </div>
        `;
    }

    renderItemDetails(item) {
        let html = '';
        
        // Define property handlers - maps property names to display titles and types
        const propertyConfig = {
            // Skip these - handled separately or not displayed
            skip: ['id', 'name', 'description', 'databases'],
            
            // Arrays - render as lists
            lists: {
                'components': 'Components',
                'proteins': 'Key Proteins',
                'lipidSynthesis': 'Lipid Synthesis Functions',
                'techniques': 'Techniques',
                'diseases': 'Associated Diseases',
                'complexes': 'Protein Complexes',
                'enzymes': 'Enzymes',
                'markers': 'Molecular Markers',
                'functions': 'Functions',
                'keyFacts': 'Key Facts',
                'features': 'Key Features',
                'substructures': 'Substructures',
                'singleCellMethods': 'Single Cell Methods',
                'nuclearBodies': 'Nuclear Bodies',
                'mechanisms': 'Mechanisms',
                'methods': 'Methods',
                'applications': 'Applications',
                'advantages': 'Advantages',
                'limitations': 'Limitations',
                'challenges': 'Challenges',
                'variants': 'Variants',
                'commercialPlatforms': 'Commercial Platforms',
                'domains': 'Domains',
                'regulation': 'Regulation',
                'biologicalInsights': 'Biological Insights',
                'biologicalRelevance': 'Biological Relevance',
                'keyAdvances': 'Key Advances',
                'singleCellAspects': 'Single Cell Aspects',
                'singleCellRelevance': 'Single Cell Relevance',
                'singleCellApproaches': 'Single Cell Approaches',
                'detection': 'Detection Methods',
                'principles': 'Principles',
                'cargo': 'Cargo Types',
                'inhibitors': 'Inhibitors',
                'imaging': 'Imaging Techniques',
                'cellularRoles': 'Cellular Roles',
                'cellularFunctions': 'Cellular Functions',
                'timescales': 'Timescales',
                'therapeuticTargets': 'Therapeutic Targets',
                'emergentProperties': 'Emergent Properties',
                'checkpoints': 'Checkpoints',
                'stages': 'Stages',
                'errors': 'Errors',
                'physiologicalRoles': 'Physiological Roles',
                'sensors': 'Sensors',
                'singleCellInsights': 'Single Cell Insights',
                'triggers': 'Triggers',
                'outcomes': 'Outcomes',
                'stressTypes': 'Stress Types',
                'analysis': 'Analysis Tools',
                'biologicalProcesses': 'Biological Processes',
                'characteristics': 'Characteristics',
                'mechanosensing': 'Mechanosensing',
                'molecularOscillators': 'Molecular Oscillators',
                'guidance': 'Guidance',
                'steps': 'Steps',
                'determinants': 'Determinants',
                'reversibility': 'Reversibility',
                'hierarchy': 'Hierarchy',
                'mechanisms': 'Mechanisms',
                'effectors': 'Effectors',
                'buffering': 'Buffering',
                'encoding': 'Encoding',
                'patterns': 'Patterns',
                'panels': 'Panels',
                'modelTypes': 'Model Types'
            },
            
            // Objects - render as generic sections
            objects: {
                'organization': 'Organization',
                'composition': 'Composition',
                'asymmetry': 'Membrane Asymmetry',
                'structure': 'Structure',
                'dynamics': 'Dynamics',
                'types': 'Types',
                'stress': 'Stress Response',
                'properties': 'Properties',
                'biogenesis': 'Biogenesis',
                'chromatin': 'Chromatin Organization',
                'keyMetrics': 'Key Metrics',
                'nuclearPoreComplex': 'Nuclear Pore Complex',
                'transport': 'Transport',
                'modifications': 'Modifications',
                'regulation': 'Regulation',
                'mechanism': 'Mechanism',
                'elongation': 'Elongation',
                'assembly': 'Assembly',
                'alternativeSplicing': 'Alternative Splicing',
                'damage': 'Damage',
                'proteins': 'Proteins',
                'lipidSynthesis': 'Lipid Synthesis',
                'marks': 'Marks',
                'structures': 'Structures',
                'cellToCell': 'Cell to Cell',
                'coverage': 'Coverage',
                'resolution': 'Resolution',
                'principle': 'Principle',
                'technology': 'Technology',
                'methodology': 'Methodology',
                'vendor': 'Vendor',
                'scale': 'Scale',
                'complexity': 'Complexity',
                'diversity': 'Diversity',
                'statistics': 'Statistics',
                'spindle': 'Spindle',
                'kinetochore': 'Kinetochore',
                'forces': 'Forces',
                'asymmetry': 'Asymmetry',
                'machinery': 'Machinery',
                'contractileRing': 'Contractile Ring',
                'positioning': 'Positioning',
                'abscission': 'Abscission',
                'phases': 'Phases',
                'singleCellVariability': 'Single Cell Variability',
                'systems': 'Systems',
                'adhesions': 'Adhesions',
                'forceSensitivity': 'Force Sensitivity',
                'stiffnessRange': 'Stiffness Range',
                'responses': 'Responses',
                'polarity': 'Polarity',
                'shapeMaintenance': 'Shape Maintenance',
                'gapJunctions': 'Gap Junctions',
                'tunnelingNanotubes': 'Tunneling Nanotubes',
                'notchSignaling': 'Notch Signaling',
                'immuneSynapse': 'Immune Synapse',
                'molecules': 'Molecules',
                'diffusion': 'Diffusion',
                'ECMinteractions': 'ECM Interactions',
                'cargo': 'Cargo',
                'uptake': 'Uptake',
                'concentration': 'Concentration',
                'vesicleRelease': 'Vesicle Release',
                'plasticity': 'Plasticity',
                'waddingtonLandscape': 'Waddington Landscape',
                'hierarchies': 'Hierarchies',
                'apoptosisRegulation': 'Apoptosis Regulation',
                'integratedStressResponse': 'Integrated Stress Response',
                'hormesis': 'Hormesis',
                'SASP': 'SASP',
                'roles': 'Roles',
                'metabolicStates': 'Metabolic States',
                'integration': 'Integration',
                'metabolicCheckpoints': 'Metabolic Checkpoints'
            }
        };
        
        // Process each property
        Object.keys(item).forEach(key => {
            // Skip if in skip list
            if (propertyConfig.skip.includes(key)) return;
            
            const value = item[key];
            if (!value) return;
            
            // Handle workflow specially
            if (key === 'workflow') {
                html += this.renderWorkflow(value);
            }
            // Handle technologies/platforms with complex nested structure
            else if ((key === 'technologies' || key === 'platforms') && typeof value === 'object' && !Array.isArray(value)) {
                html += this.renderNestedObject(propertyConfig.objects[key] || this.formatLabel(key), value);
            }
            // Handle pathways with special structure
            else if (key === 'pathways' && typeof value === 'object' && !Array.isArray(value)) {
                html += this.renderPathways(value);
            }
            // Handle arrays
            else if (Array.isArray(value)) {
                const title = propertyConfig.lists[key] || this.formatLabel(key);
                html += this.renderList(title, value, 'primary');
            }
            // Handle objects
            else if (typeof value === 'object' && !Array.isArray(value)) {
                const title = propertyConfig.objects[key] || this.formatLabel(key);
                html += this.renderGenericObject(title, value);
            }
            // Handle simple values
            else if (typeof value === 'string' || typeof value === 'number') {
                // Don't render standalone strings/numbers unless they're special
                if (key === 'concentration' || key === 'coverage' || key === 'vendor' || key === 'cellToCell') {
                    html += `<div class="info-item"><strong>${this.formatLabel(key)}:</strong> ${value}</div>`;
                }
            }
        });
        
        // Handle databases last
        if (item.databases) {
            html += this.renderDatabases(item.databases);
        }
        
        return html;
    }

    // Generic object renderer
    renderGenericObject(title, data) {
        return `
            <div class="generic-section">
                <h3>${title}</h3>
                <div class="generic-grid">
                    ${Object.entries(data).map(([key, value]) => `
                        <div class="generic-item">
                            <div class="generic-label">${this.formatLabel(key)}</div>
                            <div class="generic-value">${
                                typeof value === 'object' && !Array.isArray(value) ? 
                                    this.renderSubObject(value) : 
                                    Array.isArray(value) ? 
                                        value.join(', ') : 
                                        value
                            }</div>
                        </div>
                    `).join('')}
                </div>
            </div>
        `;
    }

    // Renderer for nested objects (like technologies, platforms)
    renderNestedObject(title, data) {
        return `
            <div class="nested-section">
                <h3>${title}</h3>
                <div class="nested-categories">
                    ${Object.entries(data).map(([category, items]) => `
                        <div class="nested-category">
                            <h4>${category}</h4>
                            ${this.renderNestedItems(items)}
                        </div>
                    `).join('')}
                </div>
            </div>
        `;
    }

    renderNestedItems(items) {
        if (Array.isArray(items)) {
            return `<ul>${items.map(item => `<li>${item}</li>`).join('')}</ul>`;
        }
        
        if (typeof items === 'object') {
            return `
                <dl class="nested-details">
                    ${Object.entries(items).map(([key, value]) => `
                        <dt>${this.formatLabel(key)}</dt>
                        <dd>${Array.isArray(value) ? value.join(', ') : value}</dd>
                    `).join('')}
                </dl>
            `;
        }
        
        return `<p>${items}</p>`;
    }

    // Handle sub-objects within generic objects
    renderSubObject(obj) {
        if (typeof obj === 'object' && !Array.isArray(obj)) {
            return Object.entries(obj).map(([k, v]) => 
                `<div><strong>${this.formatLabel(k)}:</strong> ${v}</div>`
            ).join('');
        }
        return obj;
    }

    // Special renderer for pathways (keeps complex structure)
    renderPathways(pathways) {
        return `
            <div class="pathways">
                <h3>Pathways</h3>
                ${Object.entries(pathways).map(([name, details]) => `
                    <div class="pathway-section">
                        <h4>${name}</h4>
                        ${typeof details === 'object' ?
                            `<div class="pathway-details">
                                ${Object.entries(details).map(([key, value]) => `
                                    <div class="pathway-row">
                                        <span class="pathway-key">${this.formatLabel(key)}:</span>
                                        <span class="pathway-value">${
                                            Array.isArray(value) ? value.join(', ') : value
                                        }</span>
                                    </div>
                                `).join('')}
                            </div>` :
                            `<p>${details}</p>`
                        }
                    </div>
                `).join('')}
            </div>
        `;
    }
    
    renderWorkflow(workflow) {
        return `
            <div class="workflow">
                <h3>Workflow</h3>
                <div class="workflow-steps">
                    ${workflow.map((step, index) => `
                        <div class="workflow-step">
                            <div class="step-number">${index + 1}</div>
                            <div class="step-content">${step}</div>
                        </div>
                    `).join('')}
                </div>
            </div>
        `;
    }
    
    renderList(title, items, type = 'default') {
        if (!items || (Array.isArray(items) && items.length === 0)) return '';
        
        return `
            <div class="list-section ${type}">
                <h3>${title}</h3>
                <ul>
                    ${items.map(item => `<li>${item}</li>`).join('')}
                </ul>
            </div>
        `;
    }
    
    renderDatabases(databases) {
        return `
            <div class="databases">
                <h3>Related Databases</h3>
                <div class="database-list">
                    ${databases.map(db => `
                        <a href="${db.url}" target="_blank" class="database-item">
                            <i class="fas fa-database"></i>
                            <div>
                                <div class="database-name">${db.name}</div>
                                <div class="database-description">${db.description || ''}</div>
                            </div>
                            <i class="fas fa-external-link-alt"></i>
                        </a>
                    `).join('')}
                </div>
            </div>
        `;
    }
    
    renderResources(resources) {
        if (!resources) return '';
        
        return `
            <div class="resources">
                <h2>Resources</h2>
                <div class="resource-grid">
                    ${resources.databases ? `
                        <div class="resource-section">
                            <h3>Databases</h3>
                            ${resources.databases.map(db => `
                                <a href="${db.url}" target="_blank" class="resource-link">
                                    ${db.name} - ${db.type}
                                </a>
                            `).join('')}
                        </div>
                    ` : ''}
                    
                    ${resources.tools ? `
                        <div class="resource-section">
                            <h3>Tools</h3>
                            ${resources.tools.map(tool => `
                                <a href="${tool.url}" target="_blank" class="resource-link">
                                    ${tool.name} - ${tool.purpose}
                                </a>
                            `).join('')}
                        </div>
                    ` : ''}
                    
                    ${resources.tutorials ? `
                        <div class="resource-section">
                            <h3>Tutorials</h3>
                            ${resources.tutorials.map(tutorial => `
                                <a href="${tutorial.url}" target="_blank" class="resource-link">
                                    ${tutorial.title} - ${tutorial.description}
                                </a>
                            `).join('')}
                        </div>
                    ` : ''}
                    
                    ${resources.reviews ? `
                        <div class="resource-section">
                            <h3>Key Reviews</h3>
                            ${resources.reviews.map(review => `
                                <div class="resource-link">
                                    ${review.title} (${review.year}) - ${review.journal}
                                    ${review.doi ? `<br><small>DOI: ${review.doi}</small>` : ''}
                                </div>
                            `).join('')}
                        </div>
                    ` : ''}
                </div>
            </div>
        `;
    }
    
    renderSidebar(sectionData, activeItemId) {
        return sectionData.categories.map(category => `
            <div class="sidebar-section">
                <div class="sidebar-title">${category.name}</div>
                <ul class="sidebar-list">
                    ${category.items.map(item => `
                        <li class="sidebar-item ${item.id === activeItemId ? 'active' : ''}"
                            onclick="app.navigateToSection('${sectionData.id}', '${item.id}')">
                            ${item.name}
                        </li>
                    `).join('')}
                </ul>
            </div>
        `).join('');
    }
    
    renderEmpty() {
        this.container.innerHTML = `
            <div class="empty-state">
                <i class="fas fa-inbox"></i>
                <h2>No Content Available</h2>
                <p>Please select a section from the navigation above.</p>
            </div>
        `;
    }
    
    attachCardHandlers() {
        document.querySelectorAll('.item-card').forEach(card => {
            card.addEventListener('click', () => {
                const sectionId = card.dataset.section;
                const itemId = card.dataset.item;
                this.app.navigateToSection(sectionId, itemId);
            });
        });
    }
    
    formatLabel(key) {
        return key.replace(/([A-Z])/g, ' $1')
                  .replace(/_/g, ' ')
                  .replace(/^./, str => str.toUpperCase());
    }
}