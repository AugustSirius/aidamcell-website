import { cellularComponents } from './cellular-components.js';
import { multiOmicsTechnologies } from './multi-omics-technologies.js';
import { cellularDynamics } from './cellular-dynamics.js';
import { imagingSpatial } from './imaging-spatial.js';
import { computationalTools } from './computational-tools.js';

export class KnowledgeBase {
    constructor() {
        // Step 1: Create sections with actual data first
        this.sections = {
            'cellular-components': cellularComponents,
            'multi-omics-technologies': multiOmicsTechnologies,
            'cellular-dynamics': cellularDynamics,
            'imaging-spatial': imagingSpatial,
            'computational-tools': computationalTools
        };
        
        // Step 2: Now generate overview when sections exist
        this.sections['overview'] = this.generateOverview();
    }
    
    generateOverview() {
        return {
            id: 'overview',
            name: 'SingleCellVerse - Virtual Cell Knowledge Hub',
            description: 'Comprehensive knowledge base for building AI Virtual Cell (AIVC) models through integrated single-cell biology',
            stats: {
                totalItems: this.countItems(),
                databases: this.countDatabases(),
                tools: this.countTools(),
                lastUpdated: new Date().toLocaleDateString()
            },
            highlights: [
                {
                    title: 'Cellular Components',
                    count: this.sections['cellular-components'].categories.length,
                    icon: 'fas fa-puzzle-piece',
                    color: 'primary'
                },
                {
                    title: 'Multi-Omics Technologies',
                    count: this.sections['multi-omics-technologies'].categories.length,
                    icon: 'fas fa-layer-group',
                    color: 'secondary'
                },
                {
                    title: 'Cellular Dynamics',
                    count: this.sections['cellular-dynamics'].categories.length,
                    icon: 'fas fa-sync-alt',
                    color: 'info'
                },
                {
                    title: 'Imaging & Spatial',
                    count: this.sections['imaging-spatial'].categories.length,
                    icon: 'fas fa-microscope',
                    color: 'success'
                },
                {
                    title: 'Computational Tools',
                    count: this.sections['computational-tools'].categories.length,
                    icon: 'fas fa-laptop-code',
                    color: 'warning'
                }
            ]
        };
    }

    
    getSection(sectionId) {
        return this.sections[sectionId] || null;
    }
    
    getItem(sectionId, itemId) {
        const section = this.getSection(sectionId);
        if (!section || !section.categories) return null;
        
        for (const category of section.categories) {
            const item = category.items?.find(i => i.id === itemId);
            if (item) return item;
        }
        return null;
    }
    
    search(query) {
        const results = [];
        const searchTerm = query.toLowerCase();
        
        Object.values(this.sections).forEach(section => {
            if (section.categories) {
                section.categories.forEach(category => {
                    category.items?.forEach(item => {
                        // 扩展搜索范围，包含更多字段
                        const searchableText = [
                            item.name,
                            item.description,
                            // 新增搜索字段
                            JSON.stringify(item.specifications),
                            JSON.stringify(item.features),
                            JSON.stringify(item.applications),
                            JSON.stringify(item.methods),
                            JSON.stringify(item.technologies),
                            JSON.stringify(item.platforms),
                            JSON.stringify(item.workflow),
                            JSON.stringify(item.advantages),
                            JSON.stringify(item.limitations),
                            JSON.stringify(item.keyMetrics),
                            JSON.stringify(item.keyFacts),
                            JSON.stringify(item.biologicalInsights),
                            JSON.stringify(item.challenges)
                        ].filter(Boolean).join(' ').toLowerCase();
                        
                        if (searchableText.includes(searchTerm)) {
                            results.push({
                                ...item,
                                section: section.name,
                                category: category.name,
                                sectionId: section.id
                            });
                        }
                    });
                });
            }
        });
        
        return results;
    }
    
    countItems() {
        let count = 0;
        Object.values(this.sections).forEach(section => {
            if (section.categories) {
                section.categories.forEach(category => {
                    count += category.items?.length || 0;
                });
            }
        });
        return count;
    }
    
    countDatabases() {
        const databases = new Set();
        Object.values(this.sections).forEach(section => {
            section.resources?.databases?.forEach(db => {
                databases.add(db.name);
            });
            // Also count databases in items
            section.categories?.forEach(category => {
                category.items?.forEach(item => {
                    item.databases?.forEach(db => {
                        databases.add(db.name);
                    });
                });
            });
        });
        return databases.size;
    }
    
    countTools() {
        const tools = new Set();
        Object.values(this.sections).forEach(section => {
            section.resources?.tools?.forEach(tool => {
                tools.add(tool.name);
            });
            // Also count tools/analysis software mentioned in items
            section.categories?.forEach(category => {
                category.items?.forEach(item => {
                    item.analysis?.forEach(tool => {
                        tools.add(tool);
                    });
                });
            });
        });
        return tools.size;
    }
}