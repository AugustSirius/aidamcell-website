export class Navigation {
    constructor(app) {
        this.app = app;
        this.container = document.getElementById('main-nav');
        
        this.sections = [
            {
                id: 'overview',
                name: 'Overview',
                icon: 'fas fa-home'
            },
            {
                id: 'cellular-components',
                name: 'Cellular Components',
                icon: 'fas fa-puzzle-piece'
            },
            {
                id: 'multi-omics-technologies',
                name: 'Multi-Omics Technologies',
                icon: 'fas fa-layer-group'
            },
            {
                id: 'cellular-dynamics',
                name: 'Cellular Dynamics',
                icon: 'fas fa-sync-alt'
            },
            {
                id: 'imaging-spatial',
                name: 'Imaging & Spatial',
                icon: 'fas fa-microscope'
            },
            {
                id: 'computational-tools',
                name: 'Computational Tools',
                icon: 'fas fa-laptop-code'
            }
        ];
    }
    
    render() {
        const navContainer = document.createElement('div');
        navContainer.className = 'nav-container';
        
        this.sections.forEach(section => {
            const navItem = this.createNavItem(section);
            navContainer.appendChild(navItem);
        });
        
        this.container.appendChild(navContainer);
    }
    
    createNavItem(section) {
        const item = document.createElement('div');
        item.className = 'nav-item';
        item.dataset.section = section.id;
        
        item.innerHTML = `
            <i class="${section.icon}"></i>
            <span>${section.name}</span>
        `;
        
        item.addEventListener('click', () => {
            this.app.navigateToSection(section.id);
        });
        
        return item;
    }
    
    setActive(sectionId) {
        // Remove all active classes
        document.querySelectorAll('.nav-item').forEach(item => {
            item.classList.remove('active');
        });
        
        // Add active class to current section
        const activeItem = document.querySelector(`[data-section="${sectionId}"]`);
        if (activeItem) {
            activeItem.classList.add('active');
        }
    }
}