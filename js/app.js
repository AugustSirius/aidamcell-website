// Main Application Module
import { Navigation } from './components/navigation.js';
import { Search } from './components/search.js';
import { ContentViewer } from './components/content-viewer.js';
import { KnowledgeBase } from './data/knowledge-base.js';

class SingleCellApp {
    constructor() {
        console.log('SingleCellApp initializing...');
        
        this.currentSection = 'overview';
        this.currentItem = null;
        
        try {
            this.knowledgeBase = new KnowledgeBase();
            console.log('Knowledge base loaded');
            
            this.navigation = new Navigation(this);
            this.search = new Search(this);
            this.contentViewer = new ContentViewer(this);
            
            this.init();
        } catch (error) {
            console.error('Error initializing app:', error);
        }
    }
    
    init() {
        console.log('Initializing components...');
        
        // Initialize components
        this.navigation.render();
        this.search.render();
        
        // Load initial content
        this.navigateToSection('overview');
        
        // Set up event listeners
        this.setupEventListeners();
        
        console.log('App initialized successfully');
    }
    
    setupEventListeners() {
        // Handle browser back/forward
        window.addEventListener('popstate', (e) => {
            if (e.state) {
                this.navigateToSection(e.state.section, e.state.item, false);
            }
        });
        
        // Handle hash changes (for direct URL access)
        if (window.location.hash) {
            const hash = window.location.hash.substring(1);
            const [section, item] = hash.split('/');
            if (section) {
                this.navigateToSection(section, item || null, false);
            }
        }
    }
    
    navigateToSection(sectionId, itemId = null, pushState = true) {
        console.log(`Navigating to section: ${sectionId}, item: ${itemId}`);
        
        this.currentSection = sectionId;
        this.currentItem = itemId;
        
        // Update navigation
        this.navigation.setActive(sectionId);
        
        // Load content
        const sectionData = this.knowledgeBase.getSection(sectionId);
        if (sectionData) {
            this.contentViewer.render(sectionData, itemId);
        } else {
            console.error(`Section not found: ${sectionId}`);
            this.contentViewer.renderEmpty();
        }
        
        // Update URL
        if (pushState) {
            const url = itemId ? `#${sectionId}/${itemId}` : `#${sectionId}`;
            history.pushState(
                { section: sectionId, item: itemId },
                '',
                url
            );
        }
    }
    
    search(query) {
        console.log(`Searching for: ${query}`);
        const results = this.knowledgeBase.search(query);
        console.log(`Found ${results.length} results`);
        // The search component handles display directly
        return results;
    }
}

// Initialize app when DOM is ready
document.addEventListener('DOMContentLoaded', () => {
    console.log('DOM loaded, initializing app...');
    
    try {
        window.app = new SingleCellApp();
        console.log('App created and attached to window.app');
    } catch (error) {
        console.error('Failed to initialize app:', error);
        
        // Show error message to user
        document.getElementById('content-viewer').innerHTML = `
            <div style="padding: 2rem; text-align: center; color: red;">
                <h2>Error Loading Application</h2>
                <p>Please check the browser console for details.</p>
                <p>Error: ${error.message}</p>
            </div>
        `;
    }
});

// Export for debugging
export { SingleCellApp };