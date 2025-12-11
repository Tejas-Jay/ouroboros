import requests
from bs4 import BeautifulSoup #type: ignore
import pandas as pd #type: ignore
import json
from typing import Dict, List, Optional
import time
import logging
from urllib.parse import urljoin
import re
import streamlit as st #type: ignore

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class ToxicityDatabaseScraper:
    """
    Scrapes toxicity data from multiple sources for snake venoms
    Optimized for Streamlit caching
    """
    
    def __init__(self):
        self.headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
        }
        self.session = requests.Session()
        self.session.headers.update(self.headers)
    
    def scrape_pubmed_ld50(self, snake_name: str) -> List[Dict]:
        """
        Scrape LD50 data from PubMed (via EUtils API - no API key needed)
        """
        try:
            base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            
            params = {
                'db': 'pubmed',
                'term': f'"{snake_name}" LD50 toxicity venom',
                'retmax': 10,
                'rettype': 'json'
            }
            
            response = self.session.get(base_url, params=params, timeout=10)
            response.raise_for_status()
            
            results = response.json()
            pubmed_ids = results.get('esearchresult', {}).get('idlist', [])
            
            ld50_data = []
            for pmid in pubmed_ids[:3]:
                summary = self.get_pubmed_summary(pmid)
                if summary:
                    ld50_data.append(summary)
            
            return ld50_data
            
        except Exception as e:
            logger.error(f"Error scraping PubMed for {snake_name}: {str(e)}")
            return []
    
    def get_pubmed_summary(self, pmid: str) -> Optional[Dict]:
        """
        Get article summary from PubMed
        """
        try:
            base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            
            params = {
                'db': 'pubmed',
                'id': pmid,
                'rettype': 'json'
            }
            
            response = self.session.get(base_url, params=params, timeout=10)
            response.raise_for_status()
            
            data = response.json()
            result = data.get('result', {}).get(pmid, {})
            
            return {
                'pmid': pmid,
                'title': result.get('title', 'N/A'),
                'pubdate': result.get('pubdate', 'N/A'),
                'source': 'PubMed'
            }
        except Exception as e:
            logger.error(f"Error getting PubMed summary for {pmid}: {str(e)}")
            return None
    
    def scrape_wikipedia_snake_data(self, snake_name: str) -> Dict:
        """
        Scrape snake information from Wikipedia
        """
        try:
            url = f"https://en.wikipedia.org/wiki/{snake_name.replace(' ', '_')}"
            response = self.session.get(url, timeout=10)
            response.raise_for_status()
            
            soup = BeautifulSoup(response.content, 'html.parser')
            
            data = {
                'name': snake_name,
                'description': '',
                'venom_yield': 'Data pending',
                'toxicity_level': 'Data pending',
            }
            
            # Extract infobox data
            infobox = soup.find('table', {'class': 'infobox'})
            if infobox:
                rows = infobox.find_all('tr')
                for row in rows:
                    cells = row.find_all(['td', 'th'])
                    if len(cells) >= 2:
                        key = cells[0].get_text(strip=True).lower()
                        value = cells[1].get_text(strip=True)
                        
                        if 'venom' in key:
                            data['venom_yield'] = value
                        elif 'toxicity' in key or 'lethal' in key:
                            data['toxicity_level'] = value
            
            # Extract first paragraph as description
            first_para = soup.find('p')
            if first_para:
                data['description'] = first_para.get_text(strip=True)[:300]
            
            return data
            
        except Exception as e:
            logger.error(f"Error scraping Wikipedia for {snake_name}: {str(e)}")
            return {'name': snake_name, 'description': 'Wikipedia data unavailable'}
    
    def get_clinical_effects_database(self) -> Dict:
        """
        Database of clinical effects by snake species
        """
        return {
            'Black Mamba': {
                'primary_effects': ['Neurotoxicity', 'Cardiovascular collapse'],
                'onset_time': '15-20 minutes',
                'mortality_rate': '100% if untreated',
                'treatment': 'Polyvalent antivenom (African)',
                'ld50_mice': '0.04-0.12 mg/kg'
            },
            'King Cobra': {
                'primary_effects': ['Neurotoxicity', 'Respiratory paralysis'],
                'onset_time': '30 minutes to 2 hours',
                'mortality_rate': '60-90% if untreated',
                'treatment': 'Polyvalent antivenom (Asian)',
                'ld50_mice': '0.37-0.5 mg/kg'
            },
            'Inland Taipan': {
                'primary_effects': ['Neurotoxicity', 'Coagulopathy'],
                'onset_time': '30-45 minutes',
                'mortality_rate': '100% if untreated',
                'treatment': 'Polyvalent antivenom (Australian)',
                'ld50_mice': '0.025 mg/kg'
            },
            'Indian Cobra': {
                'primary_effects': ['Neurotoxicity', 'Cardiovascular effects'],
                'onset_time': '30 minutes to 4 hours',
                'mortality_rate': '10-20% with treatment',
                'treatment': 'Polyvalent antivenom (Indian)',
                'ld50_mice': '0.8 mg/kg'
            },
            'Saw-scaled Viper': {
                'primary_effects': ['Hemorrhage', 'Coagulopathy', 'Acute kidney injury'],
                'onset_time': '1-6 hours',
                'mortality_rate': '1-25% depending on region',
                'treatment': 'Polyvalent antivenom + supportive care',
                'ld50_mice': '0.12 mg/kg'
            },
            'Puff Adder': {
                'primary_effects': ['Hemorrhage', 'Tissue necrosis'],
                'onset_time': '1-2 hours',
                'mortality_rate': '10-20% if untreated',
                'treatment': 'Polyvalent antivenom (African)',
                'ld50_mice': '0.4-0.6 mg/kg'
            },
            'Rattlesnake': {
                'primary_effects': ['Hemorrhage', 'Coagulopathy', 'Tissue damage'],
                'onset_time': '1-3 hours',
                'mortality_rate': '5-10% if untreated',
                'treatment': 'CroFab or Anavip antivenom',
                'ld50_mice': '0.5-1.2 mg/kg'
            }
        }
    
    def scrape_toxin_function_database(self) -> Dict:
        """
        Database of common snake venom components
        """
        return {
            'Serine Proteases': {
                'description': 'Break down blood proteins and disrupt coagulation',
                'clinical_effects': 'Coagulopathy, hemorrhage, tissue damage',
                'severity': 'High'
            },
            'Metalloproteinases': {
                'description': 'Degrade extracellular matrix proteins',
                'clinical_effects': 'Hemorrhage, tissue destruction, edema',
                'severity': 'High'
            },
            'Phospholipase A2': {
                'description': 'Breaks down cell membranes',
                'clinical_effects': 'Anticoagulant, hemolysis, myotoxicity, neurotoxicity',
                'severity': 'Very High'
            },
            'Three-Finger Toxins': {
                'description': 'Block neuromuscular transmission',
                'clinical_effects': 'Paralysis, respiratory failure',
                'severity': 'Very High'
            },
            'Serine Protease Inhibitors': {
                'description': 'Inhibit blood clotting factors',
                'clinical_effects': 'Bleeding disorders',
                'severity': 'High'
            },
            'Hyaluronidases': {
                'description': 'Spread toxins through tissues',
                'clinical_effects': 'Increased tissue penetration',
                'severity': 'Medium'
            }
        }
    
    def get_global_snakebite_stats(self) -> Dict:
        """
        Global snakebite statistics
        """
        return {
            'global_deaths_per_year': '81,000 - 138,000',
            'global_amputations_per_year': '400,000',
            'most_dangerous_regions': [
                'Sub-Saharan Africa',
                'South Asia',
                'Southeast Asia',
                'Latin America'
            ]
        }
    
    def compile_toxicity_report(self, snake_name: str) -> Dict:
        """
        Compile comprehensive toxicity report for a snake
        """
        clinical_db = self.get_clinical_effects_database()
        
        # Check if snake is in database
        if snake_name in clinical_db:
            clinical_data = clinical_db[snake_name]
        else:
            clinical_data = {'status': 'Data not in curated database'}
        
        report = {
            'snake_name': snake_name,
            'wikipedia_data': self.scrape_wikipedia_snake_data(snake_name),
            'pubmed_articles': self.scrape_pubmed_ld50(snake_name),
            'clinical_effects': clinical_data,
            'toxin_types': self.scrape_toxin_function_database(),
            'global_stats': self.get_global_snakebite_stats()
        }
        
        return report


@st.cache_data(ttl=86400)
def get_toxicity_report(snake_name: str) -> Dict:
    """
    Cached wrapper for toxicity report
    24 hour cache
    """
    scraper = ToxicityDatabaseScraper()
    return scraper.compile_toxicity_report(snake_name)


@st.cache_data(ttl=86400)
def get_clinical_effects() -> Dict:
    """
    Get cached clinical effects database
    """
    scraper = ToxicityDatabaseScraper()
    return scraper.get_clinical_effects_database()


@st.cache_data(ttl=86400)
def get_toxin_types() -> Dict:
    """
    Get cached toxin types database
    """
    scraper = ToxicityDatabaseScraper()
    return scraper.scrape_toxin_function_database()


@st.cache_data(ttl=86400)
def get_global_stats() -> Dict:
    """
    Get cached global snakebite statistics
    """
    scraper = ToxicityDatabaseScraper()
    return scraper.get_global_snakebite_stats()