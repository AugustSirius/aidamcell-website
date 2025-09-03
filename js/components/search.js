export class Search {
    constructor(app) {
        this.app = app;
        this.container = document.getElementById('search-container');
    }
    
    render() {
        const wrapper = document.createElement('div');
        wrapper.className = 'search-wrapper';
        wrapper.innerHTML = `
            <div class="container">
                <div class="search-box">
                    <i class="fas fa-search search-icon"></i>
                    <input 
                        type="text" 
                        id="search-input" 
                        class="search-input" 
                        placeholder="Search for components, methods, tools..."
                    >
                    <div id="search-results" class="search-results"></div>
                </div>
            </div>
        `;
        
        this.container.appendChild(wrapper);
        
        // Add event listener
        const searchInput = document.getElementById('search-input');
        searchInput.addEventListener('input', (e) => {
            this.handleSearch(e.target.value);
        });
    }
    
    handleSearch(query) {
        if (query.length < 2) {
            this.clearResults();
            return;
        }
        
        const results = this.app.knowledgeBase.search(query);
        this.displayResults(results);
    }
    
    displayResults(results) {
        const resultsContainer = document.getElementById('search-results');
        
        if (results.length === 0) {
            resultsContainer.innerHTML = '<div class="no-results">No results found</div>';
            resultsContainer.style.display = 'block';
            return;
        }
        
        resultsContainer.innerHTML = results.slice(0, 5).map(result => `
            <div class="search-result-item" data-section="${result.sectionId}" data-item="${result.id}">
                <div class="result-title">${result.name}</div>
                <div class="result-meta">${result.section} > ${result.category}</div>
            </div>
        `).join('');
        
        resultsContainer.style.display = 'block';
        
        // Add click handlers
        resultsContainer.querySelectorAll('.search-result-item').forEach(item => {
            item.addEventListener('click', () => {
                const section = item.dataset.section;
                const itemId = item.dataset.item;
                this.app.navigateToSection(section, itemId);
                this.clearResults();
                document.getElementById('search-input').value = '';
            });
        });
    }
    
    clearResults() {
        const resultsContainer = document.getElementById('search-results');
        resultsContainer.innerHTML = '';
        resultsContainer.style.display = 'none';
    }
}