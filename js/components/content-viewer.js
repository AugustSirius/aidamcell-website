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
        // 1. pick the most suitable field
        let rawTags =
            item.tags ??
            item.components ??
            item.platforms ??
            item.methods ??
            [];
    
        // 2. Normalise to an array of strings
        if (Array.isArray(rawTags)) {
            // ok
        } else if (typeof rawTags === 'object' && rawTags !== null) {
            // take the object keys (e.g. Droplet-based, Plate-based …)
            rawTags = Object.keys(rawTags);
        } else if (typeof rawTags === 'string') {
            rawTags = [rawTags];
        } else {
            rawTags = [];
        }
    
        // 3. Build the card
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
        
        // Render different properties based on what exists
        if (item.specifications) {
            html += this.renderSpecifications(item.specifications);
        }
        
        if (item.technologies) {
            html += this.renderTechnologies(item.technologies);
        }
        
        if (item.platforms) {
            html += this.renderPlatforms(item.platforms);
        }
        
        if (item.methods) {
            html += this.renderMethods(item.methods);
        }
        
        if (item.features) {
            html += this.renderList('Key Features', item.features, 'primary');
        }
        
        if (item.workflow) {
            html += this.renderWorkflow(item.workflow);
        }
        
        if (item.advantages) {
            html += this.renderList('Advantages', item.advantages, 'success');
        }
        
        if (item.limitations) {
            html += this.renderList('Limitations', item.limitations, 'warning');
        }
        
        if (item.applications) {
            html += this.renderList('Applications', item.applications, 'info');
        }
        
        if (item.biologicalInsights) {
            html += this.renderList('Biological Insights', item.biologicalInsights, 'primary');
        }
        
        if (item.challenges) {
            html += this.renderList('Challenges', item.challenges, 'warning');
        }
        
        if (item.keyMetrics) {
            html += this.renderKeyMetrics(item.keyMetrics);
        }
        
        if (item.keyFacts) {
            html += this.renderList('Key Facts', item.keyFacts, 'info');
        }
        
        if (item.keyAdvances) {
            html += this.renderList('Key Advances', item.keyAdvances, 'success');
        }
        
        if (item.databases) {
            html += this.renderDatabases(item.databases);
        }
        
        // Handle nested/complex structures
        if (item.functions) {
            html += this.renderFunctions(item.functions);
        }
        
        if (item.pathways) {
            html += this.renderPathways(item.pathways);
        }
        
        if (item.variants) {
            html += this.renderList('Variants', item.variants, 'secondary');
        }
        
        if (item.commercialPlatforms) {
            html += this.renderList('Commercial Platforms', item.commercialPlatforms, 'primary');
        }
        
        return html;
    }
    
    renderTechnologies(technologies) {
        if (Array.isArray(technologies)) {
            return this.renderList('Technologies', technologies, 'primary');
        }
        
        if (typeof technologies === 'object') {
            return `
                <div class="technologies">
                    <h3>Technologies</h3>
                    <div class="tech-categories">
                        ${Object.entries(technologies).map(([category, items]) => `
                            <div class="tech-category">
                                <h4>${category}</h4>
                                ${this.renderTechItems(items)}
                            </div>
                        `).join('')}
                    </div>
                </div>
            `;
        }
        
        return '';
    }
    
    renderTechItems(items) {
        if (Array.isArray(items)) {
            return `<ul>${items.map(item => `<li>${item}</li>`).join('')}</ul>`;
        }
        
        if (typeof items === 'object') {
            return `
                <dl class="tech-details">
                    ${Object.entries(items).map(([key, value]) => `
                        <dt>${key}</dt>
                        <dd>${Array.isArray(value) ? value.join(', ') : value}</dd>
                    `).join('')}
                </dl>
            `;
        }
        
        return `<p>${items}</p>`;
    }
    
    renderPlatforms(platforms) {
        if (typeof platforms === 'object' && !Array.isArray(platforms)) {
            return `
                <div class="platforms">
                    <h3>Platforms</h3>
                    <div class="platform-categories">
                        ${Object.entries(platforms).map(([category, items]) => `
                            <div class="platform-category">
                                <h4>${category}</h4>
                                ${typeof items === 'object' && !Array.isArray(items) ?
                                    Object.entries(items).map(([name, desc]) => 
                                        `<div class="platform-item">
                                            <strong>${name}:</strong> ${desc}
                                        </div>`
                                    ).join('') :
                                    Array.isArray(items) ?
                                        `<ul>${items.map(item => `<li>${item}</li>`).join('')}</ul>` :
                                        `<p>${items}</p>`
                                }
                            </div>
                        `).join('')}
                    </div>
                </div>
            `;
        }
        
        return this.renderList('Platforms', platforms);
    }
    
    renderMethods(methods) {
        if (Array.isArray(methods)) {
            return this.renderList('Methods', methods);
        }
        
        if (typeof methods === 'object') {
            return `
                <div class="methods">
                    <h3>Methods</h3>
                    ${Object.entries(methods).map(([name, details]) => `
                        <div class="method-section">
                            <h4>${name}</h4>
                            ${typeof details === 'object' && !Array.isArray(details) ? 
                                this.renderMethodDetails(details) : 
                                Array.isArray(details) ?
                                    `<ul>${details.map(item => `<li>${item}</li>`).join('')}</ul>` :
                                    `<p>${details}</p>`
                            }
                        </div>
                    `).join('')}
                </div>
            `;
        }
        
        return '';
    }
    
    renderMethodDetails(details) {
        return `
            <div class="method-details">
                ${Object.entries(details).map(([key, value]) => `
                    <div class="detail-row">
                        <span class="detail-key">${this.formatLabel(key)}:</span>
                        <span class="detail-value">${
                            Array.isArray(value) ? value.join(', ') : value
                        }</span>
                    </div>
                `).join('')}
            </div>
        `;
    }
    
    renderKeyMetrics(metrics) {
        return `
            <div class="key-metrics">
                <h3>Key Metrics</h3>
                <div class="metrics-grid">
                    ${Object.entries(metrics).map(([key, value]) => `
                        <div class="metric-item">
                            <div class="metric-label">${this.formatLabel(key)}</div>
                            <div class="metric-value">${value}</div>
                        </div>
                    `).join('')}
                </div>
            </div>
        `;
    }
    
    renderFunctions(functions) {
        if (Array.isArray(functions)) {
            return this.renderList('Functions', functions);
        }
        
        if (typeof functions === 'object') {
            return `
                <div class="functions">
                    <h3>Functions</h3>
                    <div class="function-grid">
                        ${Object.entries(functions).map(([key, value]) => `
                            <div class="function-item">
                                <strong>${this.formatLabel(key)}:</strong> ${value}
                            </div>
                        `).join('')}
                    </div>
                </div>
            `;
        }
        
        return '';
    }
    
    renderPathways(pathways) {
        if (typeof pathways === 'object' && !Array.isArray(pathways)) {
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
                                            <span class="pathway-key">${key}:</span>
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
        
        return this.renderList('Pathways', pathways);
    }
    
    renderSpecifications(specs) {
        return `
            <div class="specifications">
                <h3>Specifications</h3>
                <div class="spec-grid">
                    ${Object.entries(specs).map(([key, value]) => `
                        <div class="spec-item">
                            <div class="spec-label">${this.formatLabel(key)}</div>
                            <div class="spec-value">${value}</div>
                        </div>
                    `).join('')}
                </div>
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